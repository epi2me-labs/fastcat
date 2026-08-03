#include "unity.h"
#include "main.h"
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "htslib/sam.h"

#include "core/fc_alignment_summary.h"
#include "core/fc_iter.h"
#include "core/fc_read.h"
#include "core/fc_failures.h"
#include "core/fc_writer.h"

typedef struct {
    uint64_t passed_reads;
} read_count_sink_state;

static void read_count_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
    const char *fname, const sam_hdr_t *hdr)
{
    (void)read;
    (void)read_ctx;
    (void)fname;
    (void)hdr;
    read_count_sink_state *state = (read_count_sink_state *)ctx;
    state->passed_reads += 1;
}

static fc_iter_opts fq_opts_defaults(void) {
    fc_iter_opts o = {0};
    o.max_length = (size_t)-1;
    return o;
}

static const char *FASTQ_FIXTURE =
    "@r1\n"
    "ACGT\n"
    "+\n"
    "!!!!\n"
    "@r2\n"
    "ACGTAC\n"
    "+\n"
    "IIIIII\n"
    "@r3\n"
    "ACGTA\n"
    "+\n"
    "55555\n";

static const char *SAM_FIXTURE =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n"
    "r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTAC\tIIIIII\n"
    "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\t55555\n";

static const char *SAM_LOWER_QS_FIXTURE =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tqs:f:5.0\n";

static const char *SAM_MISSING_QUAL_FIXTURE =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "primary\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
    "no_qual_primary\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\t*\n"
    "secondary\t256\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t*\n";

static void make_temp_gz_fastq(char *path, size_t path_len) {
    char tmpl[] = "/tmp/fastcat_iter_fastqXXXXXX";
    int fd = mkstemp(tmpl);
    TEST_ASSERT_MESSAGE(fd != -1, "mkstemp failed for fastq fixture");
    close(fd);
    gzFile gz = gzopen(tmpl, "wb");
    TEST_ASSERT_MESSAGE(gz != NULL, "gzopen failed for fastq fixture");
    TEST_ASSERT_EQUAL_INT((int)strlen(FASTQ_FIXTURE), gzwrite(gz, FASTQ_FIXTURE, strlen(FASTQ_FIXTURE)));
    TEST_ASSERT_EQUAL_INT(0, gzclose(gz));
    char renamed[64];
    FC_SNPRINTF_CHECK(renamed, sizeof(renamed), "%s.fastq.gz", tmpl);
    TEST_ASSERT_EQUAL_INT(0, rename(tmpl, renamed));
    FC_SNPRINTF_CHECK(path, path_len, "%s", renamed);
}

static void make_temp_sam_contents(const char *sam_contents, char *path, size_t path_len) {
    char tmpl[] = "/tmp/fastcat_iter_samXXXXXX";
    int fd = mkstemp(tmpl);
    TEST_ASSERT_MESSAGE(fd != -1, "mkstemp failed for sam fixture");
    FILE *fp = fdopen(fd, "w");
    TEST_ASSERT_NOT_NULL(fp);
    TEST_ASSERT_EQUAL_INT((int)strlen(sam_contents), fprintf(fp, "%s", sam_contents));
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
    char renamed[64];
    FC_SNPRINTF_CHECK(renamed, sizeof(renamed), "%s.sam", tmpl);
    TEST_ASSERT_EQUAL_INT(0, rename(tmpl, renamed));
    FC_SNPRINTF_CHECK(path, path_len, "%s", renamed);
}

static void make_temp_sam(char *path, size_t path_len) {
    make_temp_sam_contents(SAM_FIXTURE, path, path_len);
}


// shared test driver

static void _run_iter_and_assert(
    const char *path, const fc_iter_opts *opts,
    size_t exp_nseqs, size_t exp_nbases,
    size_t exp_min_len, size_t exp_max_len,
    double exp_mean_q,
    uint64_t exp_r_record_seen,
    uint64_t exp_r_record_ok,
    uint64_t exp_r_too_short,
    uint64_t exp_r_low_quality,
    uint64_t exp_f_qual_missing
) {
    fc_iter iter = {0};
    char *paths[2] = { (char *)path, NULL };
    fc_path_list input_paths = {.paths = paths, .n_paths = 1};
    fc_file_kind expected_kind = strstr(path, ".sam") || strstr(path, ".bam")
        ? FC_FILE_BAM
        : FC_FILE_FASTQ;
    TEST_ASSERT_EQUAL_INT(0, fc_iter_init(&iter, &input_paths, 0, false, expected_kind, false));

    const char *f = NULL;
    TEST_ASSERT_TRUE(fc_iter_next_file(&iter, &f));
    TEST_ASSERT_EQUAL_INT(0, fc_iter_open_record_cursor(&iter, opts));

    fc_writer_config writer_cfg = {
        .use_barcodes = false
    };
    fc_writer *writer = fc_writer_create(&writer_cfg);
    TEST_ASSERT_NOT_NULL(writer);
    read_count_sink_state count_state = { 0 };
    fc_sink count_sink = {
        .ctx = &count_state,
        .add_read = read_count_sink_add_read
    };
    fc_writer_add_sink(writer, count_sink);
    fc_writer_begin_file(writer, f);

    size_t seen = 0;
    fc_read *r = fc_read_init();
    TEST_ASSERT_NOT_NULL(r);
    fc_iter_next_status next_rc = FC_ITER_RECORD;
    while ((next_rc = fc_iter_next_record(&iter, r)) == FC_ITER_RECORD) {
        ++seen;
        fc_writer_add_read(writer, r, f, NULL);
    }
    TEST_ASSERT_EQUAL_INT(FC_ITER_EOF, next_rc);
    fc_read_destroy(r);

    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_iter_finish_file(&iter));
    const fc_file_stats *stats = fc_iter_stats(&iter);
    TEST_ASSERT_NOT_NULL(stats);
    fc_writer_end_file(writer, f, stats);

    // some basic properties that should hold regardless of the fixture contents
    TEST_ASSERT_EQUAL_UINT(seen, stats->n_seqs);
    TEST_ASSERT_TRUE(stats->n_bases >= stats->n_seqs * stats->min_length);
    TEST_ASSERT_TRUE(stats->max_length >= stats->min_length);

    // check against expected values for the fixture
    TEST_ASSERT_EQUAL_UINT(exp_nseqs, stats->n_seqs);
    TEST_ASSERT_EQUAL_UINT(exp_nbases, stats->n_bases);
    TEST_ASSERT_EQUAL_UINT(exp_min_len, stats->min_length);
    TEST_ASSERT_EQUAL_UINT(exp_max_len, stats->max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, exp_mean_q, stats->mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats->failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_record_seen, stats->failures[R_RECORD_SEEN]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_record_ok, stats->failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_too_short, stats->failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_low_quality, stats->failures[R_LOW_QUALITY]);
    TEST_ASSERT_EQUAL_UINT64(exp_f_qual_missing, stats->failures[F_QUAL_MISSING]);
    TEST_ASSERT_EQUAL_UINT64(0, stats->failures[F_UNKNOWN_ERROR]);
    uint64_t filtered_out = stats->failures[R_TOO_LONG]
        + stats->failures[R_TOO_SHORT]
        + stats->failures[R_LOW_QUALITY]
        + stats->failures[R_DUST_MASKED];
    TEST_ASSERT_EQUAL_UINT64(exp_nseqs, count_state.passed_reads);
    TEST_ASSERT_EQUAL_UINT64(stats->failures[R_RECORD_OK], count_state.passed_reads);
    TEST_ASSERT_EQUAL_UINT64(
        stats->failures[R_RECORD_SEEN],
        stats->failures[R_RECORD_OK] + filtered_out);

    fc_writer_destroy(writer);
    fc_iter_close(&iter);
}


// run on a real fixture and check that stats are accumulated as expected

FC_TEST(iter_suite, test_fastq_iter_stats_accumulate) {
    fc_iter_opts opts = fq_opts_defaults();
    _run_iter_and_assert("test/fixtures/fastq/basic/bc1.fastq.gz", &opts,
        10, 5904, 281, 878, 12.8569, 10, 10, 0, 0, 0);
}


// check stats in a controlled case

static void _test_known_stats(const char *path) {
    _run_iter_and_assert(path, NULL, 3, 15, 4, 6, 20.0, 3, 3, 0, 0, 0);
}

FC_TEST(iter_suite, test_iter_fastq_known_stats) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_known_stats(fq_path);
    remove(fq_path);
}

FC_TEST(iter_suite, test_iter_bam_known_stats) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));
    _test_known_stats(sam_path);
    remove(sam_path);
}


// filter some things in controlled cases

static void _test_case_filters_length_and_quality(const char *path) {
    fc_iter_opts opts = fq_opts_defaults();
    opts.min_length = 5;
    opts.min_qscore = 30.0f;
    _run_iter_and_assert(path, &opts,
        1, 6, 6, 6, 40.0, 3, 1, 1, 1, 0);
}

FC_TEST(iter_suite, test_iter_fastq_filters_length_and_quality) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_case_filters_length_and_quality(fq_path);
    remove(fq_path);
}

FC_TEST(iter_suite, test_iter_bam_filters_length_and_quality) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));
    _test_case_filters_length_and_quality(sam_path);
    remove(sam_path);
}

FC_TEST(iter_suite, test_iter_bam_lowercase_qs_used_for_quality_filtering) {
    char sam_path[64];
    make_temp_sam_contents(SAM_LOWER_QS_FIXTURE, sam_path, sizeof(sam_path));
    fc_iter_opts opts = fq_opts_defaults();
    opts.min_qscore = 30.0f;
    _run_iter_and_assert(sam_path, &opts, 0, 0, 0, 0, 0.0, 1, 0, 0, 1, 0);
    remove(sam_path);
}

FC_TEST(iter_suite, test_iter_bam_missing_quality_reaches_record_level_accounting_only) {
    char sam_path[64];
    make_temp_sam_contents(SAM_MISSING_QUAL_FIXTURE, sam_path, sizeof(sam_path));
    fc_iter_opts opts = fq_opts_defaults();

    fc_iter iter = { 0 };
    char *paths[] = { sam_path };
    fc_path_list input_paths = { .paths = paths, .n_paths = 1 };
    TEST_ASSERT_EQUAL_INT(0, fc_iter_init(&iter, &input_paths, 0, false, FC_FILE_BAM, false));
    const char *fname = NULL;
    TEST_ASSERT_TRUE(fc_iter_next_file(&iter, &fname));
    TEST_ASSERT_EQUAL_INT(0, fc_iter_open_record_cursor(&iter, &opts));

    fc_alignment_summary *alignment = fc_alignment_summary_create(1, false, false, 1);
    TEST_ASSERT_NOT_NULL(alignment);
    fc_writer_config writer_cfg = { .alignment_summary = alignment };
    fc_writer *writer = fc_writer_create(&writer_cfg);
    TEST_ASSERT_NOT_NULL(writer);
    read_count_sink_state count_state = { 0 };
    fc_writer_add_sink(writer, (fc_sink) {
        .ctx = &count_state,
        .add_read = read_count_sink_add_read
    });

    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    size_t records = 0;
    while (fc_iter_next_record(&iter, read) == FC_ITER_RECORD) {
        records++;
        fc_writer_add_read(writer, read, fname, iter.bmi.hdr);
    }

    const fc_file_stats *stats = fc_iter_stats(&iter);
    TEST_ASSERT_EQUAL_UINT64(3, records);
    TEST_ASSERT_EQUAL_UINT64(1, stats->n_seqs);
    TEST_ASSERT_EQUAL_UINT64(1, stats->failures[F_QUAL_MISSING]);
    TEST_ASSERT_EQUAL_UINT64(1, count_state.passed_reads);

    fc_alignment_totals totals = { 0 };
    fc_alignment_counts_materialize(fc_alignment_summary_global(alignment), &totals);
    TEST_ASSERT_EQUAL_UINT64(3, totals.record_counts[FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(2, totals.record_counts[FC_FLAGSTAT_PRIMARY]);
    TEST_ASSERT_EQUAL_UINT64(1, totals.record_counts[FC_FLAGSTAT_SECONDARY]);

    fc_read_destroy(read);
    fc_writer_destroy(writer);
    fc_alignment_summary_destroy(alignment);
    fc_iter_close(&iter);
    remove(sam_path);
}


// all filtered sets some things

static void _test_case_all_filtered(const char *path) {
    fc_iter_opts opts = fq_opts_defaults();
    opts.min_length = 100;
    _run_iter_and_assert(path, &opts, 0, 0, 0, 0, 0.0, 3, 0, 3, 0, 0);
}

FC_TEST(iter_suite, test_iter_fastq_all_filtered_sets_minmax_zero) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_case_all_filtered(fq_path);
    remove(fq_path);
}

FC_TEST(iter_suite, test_iter_bam_all_filtered_sets_minmax_zero) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));
    _test_case_all_filtered(sam_path);
    remove(sam_path);
}


// check that open failure is recorded as unknown error via unified iterator

static void _test_case_open_failure_sets_unknown_error(const char *path, fc_file_kind kind) {
    fc_iter_opts opts = fq_opts_defaults();
    fc_iter iter = {0};
    iter.expected_kind = kind;
    iter.current = path;
    int rc = fc_iter_open_record_cursor(&iter, &opts);
    TEST_ASSERT_NOT_EQUAL(0, rc);
    const fc_file_stats *stats = fc_iter_stats(&iter);
    TEST_ASSERT_NOT_NULL(stats);
    TEST_ASSERT_EQUAL_UINT64(1, stats->failures[F_UNKNOWN_ERROR]);
    fc_iter_close(&iter);
}

FC_TEST(iter_suite, test_iter_fastq_iter_open_failure_sets_unknown_error) {
    _test_case_open_failure_sets_unknown_error("/no/such/file.fastq", FC_FILE_FASTQ);
}

FC_TEST(iter_suite, test_iter_bam_iter_open_failure_sets_unknown_error) {
    _test_case_open_failure_sets_unknown_error("/no/such/file.bam", FC_FILE_BAM);
}

// go go gadget
