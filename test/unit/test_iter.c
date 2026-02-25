#include "unity.h"
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "htslib/sam.h"

#include "fc_iter.h"
#include "fc_read.h"
#include "fc_failures.h"

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
    snprintf(renamed, sizeof(renamed), "%s.fastq.gz", tmpl);
    TEST_ASSERT_EQUAL_INT(0, rename(tmpl, renamed));
    snprintf(path, path_len, "%s", renamed);
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
    snprintf(renamed, sizeof(renamed), "%s.sam", tmpl);
    TEST_ASSERT_EQUAL_INT(0, rename(tmpl, renamed));
    snprintf(path, path_len, "%s", renamed);
}

static void make_temp_sam(char *path, size_t path_len) {
    make_temp_sam_contents(SAM_FIXTURE, path, path_len);
}


// shared test driver

static void _run_iter_and_assert(
    const char *path, const fc_iter_opts *opts,
    size_t exp_nseqs, size_t exp_nbases,
    size_t exp_min_len, size_t exp_max_len,
    double exp_mean_q, uint64_t exp_r_record_ok, uint64_t exp_r_too_short, uint64_t exp_r_low_quality
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

    size_t seen = 0;
    fc_read *r = fc_read_init();
    TEST_ASSERT_NOT_NULL(r);
    fc_iter_next_status next_rc = FC_ITER_RECORD;
    while ((next_rc = fc_iter_next_record(&iter, r)) == FC_ITER_RECORD) {
        ++seen;
    }
    TEST_ASSERT_EQUAL_INT(FC_ITER_EOF, next_rc);
    fc_read_destroy(r);

    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_iter_finish_file(&iter));
    const fc_file_stats *stats = fc_iter_stats(&iter);
    TEST_ASSERT_NOT_NULL(stats);

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
    TEST_ASSERT_EQUAL_UINT64(exp_r_record_ok, stats->failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_too_short, stats->failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(exp_r_low_quality, stats->failures[R_LOW_QUALITY]);
    TEST_ASSERT_EQUAL_UINT64(0, stats->failures[F_UNKNOWN_ERROR]);

    fc_iter_close(&iter);
}


// run on a real fixture and check that stats are accumulated as expected

static void test_fastq_iter_stats_accumulate(void) {
    fc_iter_opts opts = fq_opts_defaults();
    _run_iter_and_assert("test/data/bc1.fastq.gz", &opts,
        10, 5904, 281, 878,12.8569, 10, 0, 0);
}


// check stats in a controlled case

static void _test_known_stats(const char *path) {
    _run_iter_and_assert(path, NULL, 3, 15, 4, 6, 20.0, 3, 0, 0);
}

static void test_iter_fastq_known_stats(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_known_stats(fq_path);
    remove(fq_path);
}

static void test_iter_bam_known_stats(void) {
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
        1, 6, 6, 6, 40.0, 3, 1, 1);
}

static void test_iter_fastq_filters_length_and_quality(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_case_filters_length_and_quality(fq_path);
    remove(fq_path);
}

static void test_iter_bam_filters_length_and_quality(void) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));
    _test_case_filters_length_and_quality(sam_path);
    remove(sam_path);
}

static void test_iter_bam_lowercase_qs_used_for_quality_filtering(void) {
    char sam_path[64];
    make_temp_sam_contents(SAM_LOWER_QS_FIXTURE, sam_path, sizeof(sam_path));
    fc_iter_opts opts = fq_opts_defaults();
    opts.min_qscore = 30.0f;
    _run_iter_and_assert(sam_path, &opts, 0, 0, 0, 0, 0.0, 1, 0, 1);
    remove(sam_path);
}


// all filtered sets some things

static void _test_case_all_filtered(const char *path) {
    fc_iter_opts opts = fq_opts_defaults();
    opts.min_length = 100;
    _run_iter_and_assert(path, &opts, 0, 0, 0, 0, 0.0, 3, 3, 0);
}

static void test_iter_fastq_all_filtered_sets_minmax_zero(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));
    _test_case_all_filtered(fq_path);
    remove(fq_path);
}

static void test_iter_bam_all_filtered_sets_minmax_zero(void) {
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

static void test_iter_fastq_iter_open_failure_sets_unknown_error(void) {
    _test_case_open_failure_sets_unknown_error("/no/such/file.fastq", FC_FILE_FASTQ);
}

static void test_iter_bam_iter_open_failure_sets_unknown_error(void) {
    _test_case_open_failure_sets_unknown_error("/no/such/file.bam", FC_FILE_BAM);
}

// go go gadget

void iter_suite(void) {
    // real data
    RUN_TEST(test_fastq_iter_stats_accumulate);
    // synthetic data with known stats
    RUN_TEST(test_iter_fastq_known_stats);
    RUN_TEST(test_iter_bam_known_stats);
    // filters in controlled cases
    RUN_TEST(test_iter_fastq_filters_length_and_quality);
    RUN_TEST(test_iter_bam_filters_length_and_quality);
    RUN_TEST(test_iter_bam_lowercase_qs_used_for_quality_filtering);
    // all filtered sets min/max length to zero
    RUN_TEST(test_iter_fastq_all_filtered_sets_minmax_zero);
    RUN_TEST(test_iter_bam_all_filtered_sets_minmax_zero);
    // open failure sets unknown error
    RUN_TEST(test_iter_fastq_iter_open_failure_sets_unknown_error);
    RUN_TEST(test_iter_bam_iter_open_failure_sets_unknown_error);
}
