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
#include "fastcat/parsing.h"

static fc_fastq_iter_opts fq_opts_defaults(void) {
    fc_fastq_iter_opts o = {0};
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

static void make_temp_gz_fastq(char *path, size_t path_len) {
    char tmpl[] = "/tmp/fastcat_iter_fastqXXXXXX";
    int fd = mkstemp(tmpl);
    TEST_ASSERT_MESSAGE(fd != -1, "mkstemp failed for fastq fixture");
    close(fd);
    gzFile gz = gzopen(tmpl, "wb");
    TEST_ASSERT_MESSAGE(gz != NULL, "gzopen failed for fastq fixture");
    TEST_ASSERT_EQUAL_INT((int)strlen(FASTQ_FIXTURE), gzwrite(gz, FASTQ_FIXTURE, strlen(FASTQ_FIXTURE)));
    TEST_ASSERT_EQUAL_INT(0, gzclose(gz));
    snprintf(path, path_len, "%s", tmpl);
}

static void make_temp_sam(char *path, size_t path_len) {
    char tmpl[] = "/tmp/fastcat_iter_samXXXXXX";
    int fd = mkstemp(tmpl);
    TEST_ASSERT_MESSAGE(fd != -1, "mkstemp failed for sam fixture");
    FILE *fp = fdopen(fd, "w");
    TEST_ASSERT_NOT_NULL(fp);
    TEST_ASSERT_EQUAL_INT((int)strlen(SAM_FIXTURE), fprintf(fp, "%s", SAM_FIXTURE));
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
    snprintf(path, path_len, "%s", tmpl);
}

static void test_fastq_iter_stats_accumulate(void) {
    fc_fastq_iter iter = {0};
    fc_fastq_iter_opts opts = fq_opts_defaults();
    int rc = fc_fastq_iter_open(&iter, "test/data/bc1.fastq.gz", &opts);
    TEST_ASSERT_EQUAL_INT(0, rc);

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_fastq_iter_next(&iter)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_GREATER_THAN_UINT(0, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_fastq_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(seen, stats.n_seqs);
    TEST_ASSERT_TRUE(stats.n_bases >= stats.n_seqs * stats.min_length);
    TEST_ASSERT_TRUE(stats.max_length >= stats.min_length);
}

static void test_fastq_iter_known_stats(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));

    fc_fastq_iter iter = {0};
    fc_fastq_iter_opts opts = fq_opts_defaults();
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_iter_open(&iter, fq_path, &opts));

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_fastq_iter_next(&iter)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_EQUAL_UINT(3, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_fastq_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(3, stats.n_seqs);
    TEST_ASSERT_EQUAL_UINT(15, stats.n_bases);
    TEST_ASSERT_EQUAL_UINT(4, stats.min_length);
    TEST_ASSERT_EQUAL_UINT(6, stats.max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 20.0, stats.mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(0, stats.failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(0, stats.failures[R_LOW_QUALITY]);

    remove(fq_path);
}

static void test_fastq_iter_filters_length_and_quality(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));

    fc_fastq_iter iter = {0};
    fc_fastq_iter_opts opts = fq_opts_defaults();
    opts.min_length = 5;
    opts.min_qscore = 30.0f;

    TEST_ASSERT_EQUAL_INT(0, fc_fastq_iter_open(&iter, fq_path, &opts));

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_fastq_iter_next(&iter)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_EQUAL_UINT(1, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_fastq_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(1, stats.n_seqs);
    TEST_ASSERT_EQUAL_UINT(6, stats.n_bases);
    TEST_ASSERT_EQUAL_UINT(6, stats.min_length);
    TEST_ASSERT_EQUAL_UINT(6, stats.max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 40.0, stats.mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[R_LOW_QUALITY]);

    remove(fq_path);
}

static void test_fastq_iter_open_failure_sets_unknown_error(void) {
    fc_fastq_iter iter = {0};
    fc_fastq_iter_opts opts = fq_opts_defaults();

    int rc = fc_fastq_iter_open(&iter, "/no/such/file.fastq", &opts);
    TEST_ASSERT_NOT_EQUAL(0, rc);
    TEST_ASSERT_EQUAL_UINT64(1, iter.stats.failures[F_UNKNOWN_ERROR]);
}

static void test_fastq_iter_all_filtered_sets_minmax_zero(void) {
    char fq_path[64];
    make_temp_gz_fastq(fq_path, sizeof(fq_path));

    fc_fastq_iter iter = {0};
    fc_fastq_iter_opts opts = fq_opts_defaults();
    opts.min_length = 100; // filter everything

    TEST_ASSERT_EQUAL_INT(0, fc_fastq_iter_open(&iter, fq_path, &opts));

    fc_read *r = NULL;
    while ((r = fc_fastq_iter_next(&iter)) != NULL) {
        fc_read_destroy(r);
    }

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_fastq_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(0, stats.n_seqs);
    TEST_ASSERT_EQUAL_UINT(0, stats.n_bases);
    TEST_ASSERT_EQUAL_UINT(0, stats.min_length);
    TEST_ASSERT_EQUAL_UINT(0, stats.max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.0, stats.mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_TOO_SHORT]);

    remove(fq_path);
}

static void test_bam_iter_stats_accumulate(void) {
    htsFile *fp = hts_open("test/bamstats/400ecoli.bam", "rb");
    TEST_ASSERT_NOT_NULL(fp);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_iter iter = {0};
    fc_bam_iter_opts opts = {0};
    int rc = fc_bam_iter_open(&iter, fp, NULL, hdr, NULL, 0, INT64_MAX, &opts);
    TEST_ASSERT_EQUAL_INT(0, rc);

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_bam_iter_next(&iter, hdr)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_GREATER_THAN_UINT(0, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_bam_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(seen, stats.n_seqs);
    TEST_ASSERT_TRUE(stats.max_length >= stats.min_length);

    fc_bam_iter_close(&iter);
    sam_hdr_destroy(hdr);
    hts_close(fp);
}

static void test_bam_iter_known_stats(void) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));

    htsFile *fp = hts_open(sam_path, "r");
    TEST_ASSERT_NOT_NULL(fp);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_iter iter = {0};
    fc_bam_iter_opts opts = {0};
    TEST_ASSERT_EQUAL_INT(0, fc_bam_iter_open(&iter, fp, NULL, hdr, NULL, 0, INT64_MAX, &opts));

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_bam_iter_next(&iter, hdr)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_EQUAL_UINT(3, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_bam_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(3, stats.n_seqs);
    TEST_ASSERT_EQUAL_UINT(15, stats.n_bases);
    TEST_ASSERT_EQUAL_UINT(4, stats.min_length);
    TEST_ASSERT_EQUAL_UINT(6, stats.max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 20.0, stats.mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(0, stats.failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(0, stats.failures[R_LOW_QUALITY]);

    fc_bam_iter_close(&iter);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    remove(sam_path);
}

static void test_bam_iter_filters_length_and_quality(void) {
    char sam_path[64];
    make_temp_sam(sam_path, sizeof(sam_path));

    htsFile *fp = hts_open(sam_path, "r");
    TEST_ASSERT_NOT_NULL(fp);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_iter iter = {0};
    fc_bam_iter_opts opts = {0};
    opts.min_length = 5;
    opts.min_qscore = 30.0f;
    TEST_ASSERT_EQUAL_INT(0, fc_bam_iter_open(&iter, fp, NULL, hdr, NULL, 0, INT64_MAX, &opts));

    size_t seen = 0;
    fc_read *r = NULL;
    while ((r = fc_bam_iter_next(&iter, hdr)) != NULL) {
        ++seen;
        fc_read_destroy(r);
    }
    TEST_ASSERT_EQUAL_UINT(1, seen);

    fc_file_stats stats = {0};
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, fc_bam_iter_finish(&iter, &stats));

    TEST_ASSERT_EQUAL_UINT(1, stats.n_seqs);
    TEST_ASSERT_EQUAL_UINT(6, stats.n_bases);
    TEST_ASSERT_EQUAL_UINT(6, stats.min_length);
    TEST_ASSERT_EQUAL_UINT(6, stats.max_length);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 40.0, stats.mean_quality);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[F_FILE_OK]);
    TEST_ASSERT_EQUAL_UINT64(3, stats.failures[R_RECORD_OK]);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[R_TOO_SHORT]);
    TEST_ASSERT_EQUAL_UINT64(1, stats.failures[R_LOW_QUALITY]);

    fc_bam_iter_close(&iter);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    remove(sam_path);
}

void iter_suite(void) {
    RUN_TEST(test_fastq_iter_stats_accumulate);
    RUN_TEST(test_fastq_iter_known_stats);
    RUN_TEST(test_fastq_iter_filters_length_and_quality);
    RUN_TEST(test_fastq_iter_open_failure_sets_unknown_error);
    RUN_TEST(test_fastq_iter_all_filtered_sets_minmax_zero);
    RUN_TEST(test_bam_iter_stats_accumulate);
    RUN_TEST(test_bam_iter_known_stats);
    RUN_TEST(test_bam_iter_filters_length_and_quality);
}
