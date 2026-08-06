#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "unity.h"
#include "main.h"
#include "htslib/sam.h"

#include "core/fc_coverage.h"
#include "core/fc_common.h"
#include "test_fs_utils.h"

static bam_hdr_t *make_header(const int *lengths, size_t n_lengths)
{
    char text[2048];
    size_t used = 0;

    FC_SNPRINTF_CHECK(text + used, sizeof(text) - used, "@HD\tVN:1.6\tSO:coordinate\n");
    size_t n = strlen(text + used);
    TEST_ASSERT_TRUE(n > 0);
    TEST_ASSERT_TRUE(n < sizeof(text) - used);
    used += n;

    for (size_t i = 0; i < n_lengths; ++i) {
        TEST_ASSERT_TRUE(lengths[i] > 0);
        FC_SNPRINTF_CHECK(
            text + used, sizeof(text) - used,
            "@SQ\tSN:chr%zu\tLN:%d\n", i + 1, lengths[i]);
        n = strlen(text + used);
        TEST_ASSERT_TRUE(n > 0);
        TEST_ASSERT_TRUE(n < sizeof(text) - used);
        used += n;
    }

    bam_hdr_t *hdr = sam_hdr_parse((int)used, text);
    TEST_ASSERT_NOT_NULL(hdr);
    return hdr;
}

static void cleanup_output_dir(const char *dir_path)
{
    char path[512];

    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/global.summary.txt", dir_path);
    TEST_ASSERT_EQUAL_INT(0, remove(path));

    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/global.dist.txt", dir_path);
    TEST_ASSERT_EQUAL_INT(0, remove(path));

    TEST_ASSERT_EQUAL_INT(0, rmdir(dir_path));
}

static bam1_t *make_read(
    const char *qname, int32_t tid, int32_t pos,
    uint32_t read_len, const uint32_t *cigar, uint32_t n_cigar, const uint16_t *flags)
{
    bam1_t *b = bam_init1();
    TEST_ASSERT_NOT_NULL(b);

    uint32_t fallback_cigar = 0;
    const uint32_t *use_cigar = cigar;
    uint32_t use_n_cigar = n_cigar;
    size_t use_read_len = read_len;

    if (use_cigar == NULL) {
        // tell us where you're from
        fallback_cigar = bam_cigar_gen(read_len, BAM_CMATCH);
        use_cigar = &fallback_cigar;
        use_n_cigar = 1;
    } else {
        // what you want to become
        TEST_ASSERT_TRUE(use_n_cigar > 0);
        use_read_len = 0;
        for (uint32_t i = 0; i < use_n_cigar; ++i) {
            uint32_t op = bam_cigar_op(use_cigar[i]);
            if (bam_cigar_type(op) & 1) {
                use_read_len += bam_cigar_oplen(use_cigar[i]);
            }
        }
    }
    // and we'll say if you're ok
    TEST_ASSERT_TRUE(use_read_len > 0);

    uint16_t use_flags = flags ? *flags : 0;

    char *seq = malloc(use_read_len + 1);
    char *qual = malloc(use_read_len + 1);
    TEST_ASSERT_NOT_NULL(seq);
    TEST_ASSERT_NOT_NULL(qual);
    memset(seq, 'N', use_read_len);
    memset(qual, 'A', use_read_len);
    seq[use_read_len] = '\0';
    qual[use_read_len] = '\0';

    int rc = bam_set1(
        b,
        strlen(qname) + 1, qname,
        use_flags, tid, pos, 60,
        use_n_cigar, use_cigar,
        -1, -1, 0,
        use_read_len, seq, qual,
        0);
    free(seq);
    free(qual);

    TEST_ASSERT_TRUE(rc >= 0);
    return b;
}

FC_TEST(bamcoverage_suite, test_coverage_process_exclude_flags)
{
    const int lengths[1] = { 100 };
    const size_t active_len = (size_t)lengths[0] + 1u;

    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_flags_exclXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    // setup writer with an explicit exclude flag
    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        BAM_FQCFAIL, 0, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    // process one allowed read and one excluded read
    uint16_t flags_qcfail = BAM_FQCFAIL;
    bam1_t *read_ok = make_read("ok", 0, 10, 10, NULL, 0, NULL);
    bam1_t *read_qcfail = make_read("qcfail", 0, 30, 10, NULL, 0, &flags_qcfail);
    coverage_process(w, read_ok);
    coverage_process(w, read_qcfail);

    // only the allowed read should contribute to forward diffs
    int32_t *expected = calloc(active_len, sizeof(*expected));
    TEST_ASSERT_NOT_NULL(expected);
    expected[10] = 1;
    expected[20] = -1;
    TEST_ASSERT_EQUAL_INT32_ARRAY(expected, w->diff_fwd, active_len);
    TEST_ASSERT_EACH_EQUAL_INT32(0, w->diff_rev, active_len);
    free(expected);

    bam_destroy1(read_ok);
    bam_destroy1(read_qcfail);
    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_process_include_flags)
{
    const int lengths[1] = { 100 };
    const size_t active_len = (size_t)lengths[0] + 1u;

    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_flags_inclXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    // setup writer requiring reads to match include flag
    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        0, BAM_FREVERSE, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    // process one non-matching forward read and one matching reverse read
    uint16_t flags_rev = BAM_FREVERSE;
    bam1_t *read_fwd = make_read("fwd", 0, 40, 10, NULL, 0, NULL);
    bam1_t *read_rev = make_read("rev", 0, 50, 10, NULL, 0, &flags_rev);
    coverage_process(w, read_fwd);
    coverage_process(w, read_rev);

    // only reverse buffer should be populated by the included read
    int32_t *expected_fwd = calloc(active_len, sizeof(*expected_fwd));
    int32_t *expected_rev = calloc(active_len, sizeof(*expected_rev));
    TEST_ASSERT_NOT_NULL(expected_fwd);
    TEST_ASSERT_NOT_NULL(expected_rev);
    expected_rev[50] = 1;
    expected_rev[60] = -1;
    TEST_ASSERT_EQUAL_INT32_ARRAY(expected_fwd, w->diff_fwd, active_len);
    TEST_ASSERT_EQUAL_INT32_ARRAY(expected_rev, w->diff_rev, active_len);
    free(expected_fwd);
    free(expected_rev);

    bam_destroy1(read_fwd);
    bam_destroy1(read_rev);
    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_process_use_cigar_mixed_ops)
{
    const int lengths[1] = { 100 };
    const size_t active_len = (size_t)lengths[0] + 1u;

    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_cigarXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    // setup writer in CIGAR-aware mode
    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, true,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    // 4M2I3M1D2M2N1M at POS=10. Covered ranges:
    // [10,14), [14,17), [18,20), [22,23)
    uint32_t cigar[] = {
        bam_cigar_gen(4, BAM_CMATCH),
        bam_cigar_gen(2, BAM_CINS),
        bam_cigar_gen(3, BAM_CMATCH),
        bam_cigar_gen(1, BAM_CDEL),
        bam_cigar_gen(2, BAM_CMATCH),
        bam_cigar_gen(2, BAM_CREF_SKIP),
        bam_cigar_gen(1, BAM_CMATCH),
    };
    bam1_t *read = make_read("mix", 0, 10, 0, cigar, 7, NULL);
    coverage_process(w, read);

    // expected diff is the sum of the CIGAR-derived covered segments
    int32_t *expected = calloc(active_len, sizeof(*expected));
    TEST_ASSERT_NOT_NULL(expected);
    expected[10] += 1;
    expected[17] -= 1;
    expected[18] += 1;
    expected[20] -= 1;
    expected[22] += 1;
    expected[23] -= 1;
    TEST_ASSERT_EQUAL_INT32_ARRAY(expected, w->diff_fwd, active_len);
    TEST_ASSERT_EACH_EQUAL_INT32(0, w->diff_rev, active_len);
    free(expected);

    bam_destroy1(read);
    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_process_end_marker_at_contig_len)
{
    const int lengths[1] = { 50 };
    const size_t active_len = (size_t)lengths[0] + 1u;

    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_boundaryXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    // setup a simple writer and process one boundary read
    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    bam1_t *read = make_read("boundary", 0, 40, 10, NULL, 0, NULL);
    coverage_process(w, read);

    // end marker should land exactly at contig_len
    int32_t *expected = calloc(active_len, sizeof(*expected));
    TEST_ASSERT_NOT_NULL(expected);
    expected[40] = 1;
    expected[50] = -1;
    TEST_ASSERT_EQUAL_INT32_ARRAY(expected, w->diff_fwd, active_len);
    TEST_ASSERT_EACH_EQUAL_INT32(0, w->diff_rev, active_len);
    free(expected);

    bam_destroy1(read);
    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_writer_init_destroy_no_data)
{
    const int lengths[1] = { 80 };
    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_lifecycleXXXXXX";
    char summary[512];
    char dist[512];

    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    // init writer without processing any reads
    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    // writer starts in an uninitialized contig state
    TEST_ASSERT_EQUAL_INT(-1, w->tid);
    TEST_ASSERT_EQUAL_INT64(0, w->contig_len);
    TEST_ASSERT_NULL(w->chrom);
    TEST_ASSERT_EQUAL_UINT(0u, w->buf_size);
    TEST_ASSERT_NULL(w->diff_fwd);
    TEST_ASSERT_NULL(w->diff_rev);

    // destroy should still emit expected summary/dist outputs
    destroy_coverage_writer(w);

    FC_SNPRINTF_CHECK(summary, sizeof(summary), "%s/global.summary.txt", out_dir);
    FC_SNPRINTF_CHECK(dist, sizeof(dist), "%s/global.dist.txt", out_dir);
    TEST_ASSERT_EQUAL_INT(0, access(summary, F_OK));
    TEST_ASSERT_EQUAL_INT(0, access(dist, F_OK));

    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_thresholds_deduplicated)
{
    const int lengths[1] = { 80 };
    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_thresholdsXXXXXX";
    uint32_t thresholds[] = { 20, 5, 20, 10, 5 };
    uint32_t expected[] = { 5, 10, 20 };
 
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    bam_hdr_t *hdr = make_header(lengths, 1);
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0,
        thresholds, sizeof(thresholds) / sizeof(thresholds[0]),
        NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);
    TEST_ASSERT_EQUAL_UINT(3u, w->writers[0]->n_thresholds);
    TEST_ASSERT_EQUAL_UINT32_ARRAY(expected, w->writers[0]->thresholds, 3);

    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_summary_thresholds_handle_depth_above_initial_cap)
{
    const int lengths[1] = { 100 };
    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_thresholdsXXXXXX";
    char summary_path[512];

    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    bam_hdr_t *hdr = make_header(lengths, 1);
    uint32_t thresholds[] = { 1, 5, 10, 20 };
    cov_writer w = init_coverage_writer(
        out_dir, false, true, false,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0, thresholds, 4, NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    for (int i = 0; i < 300; ++i) {
        char qname[32];
        FC_SNPRINTF_CHECK(qname, sizeof(qname), "depth_%d", i);
        bam1_t *read = make_read(qname, 0, 0, 1, NULL, 0, NULL);
        coverage_process(w, read);
        bam_destroy1(read);
    }

    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);

    FC_SNPRINTF_CHECK(summary_path, sizeof(summary_path), "%s/global.summary.txt", out_dir);
    FILE *fp = fopen(summary_path, "r");
    TEST_ASSERT_NOT_NULL(fp);

    char line[512];
    TEST_ASSERT_NOT_NULL(fgets(line, sizeof(line), fp)); // header
    TEST_ASSERT_NOT_NULL(fgets(line, sizeof(line), fp)); // chr1 row

    char chrom[64] = { 0 };
    long long start = 0, end = 0, length = 0, bases = 0, min_cov = 0, max_cov = 0;
    double mean = 0.0, t1 = 0.0, t5 = 0.0, t10 = 0.0, t20 = 0.0;
    int parsed = sscanf(line,
        "%63[^\t]\t%lld\t%lld\t%lld\t%lld\t%lf\t%lld\t%lld\t%lf\t%lf\t%lf\t%lf",
        chrom, &start, &end, &length, &bases, &mean, &min_cov, &max_cov, &t1, &t5, &t10, &t20);
    TEST_ASSERT_EQUAL_INT(12, parsed);
    TEST_ASSERT_EQUAL_STRING("chr1", chrom);
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.01f, (float)t1);
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.01f, (float)t5);
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.01f, (float)t10);
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.01f, (float)t20);

    fclose(fp);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_summary_snapshot_preserves_numeric_precision)
{
    const int lengths[1] = { 7 };
    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_snapshotXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    bam_hdr_t *hdr = make_header(lengths, 1);
    fc_coverage_summaries *summaries = fc_coverage_summaries_create(1);
    cov_writer writer = init_coverage_writer(
        out_dir, false, true, false,
        -1, -1, false,
        hdr, NULL, NULL, NULL, 0, NULL, 0, NULL, 0,
        fc_coverage_summaries_prepare(summaries, 0));
    bam1_t *read = make_read("read", 0, 0, 3, NULL, 0, NULL);
    coverage_process(writer, read);
    bam_destroy1(read);
    destroy_coverage_writer(writer);

    const fc_coverage_summary *summary = fc_coverage_summaries_get(summaries, 0);
    TEST_ASSERT_NOT_NULL(summary);
    TEST_ASSERT_EQUAL_size_t(1, summary->n_targets);
    const fc_coverage_summary_target *target = &summary->targets[0];
    TEST_ASSERT_EQUAL_STRING("global", target->name);
    TEST_ASSERT_EQUAL_size_t(2, target->n_summary_rows);
    TEST_ASSERT_EQUAL_STRING("chr1", target->summary_rows[0].chrom);
    TEST_ASSERT_EQUAL_size_t(7, target->summary_rows[0].length);
    TEST_ASSERT_EQUAL_size_t(3, target->summary_rows[0].bases);
    TEST_ASSERT_EQUAL_DOUBLE(3.0 / 7.0, target->summary_rows[0].mean);
    TEST_ASSERT_EQUAL_DOUBLE(3.0 / 7.0, target->summary_rows[0].threshold_fractions[0]);
    TEST_ASSERT_EQUAL_STRING("global", target->summary_rows[1].chrom);
    TEST_ASSERT_EQUAL_size_t(2, target->n_distribution_rows);
    TEST_ASSERT_EQUAL_size_t(1, target->distribution_rows[0].depth);
    TEST_ASSERT_EQUAL_DOUBLE(3.0 / 7.0, target->distribution_rows[0].fraction);
    TEST_ASSERT_EQUAL_size_t(0, target->distribution_rows[1].depth);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, target->distribution_rows[1].fraction);

    fc_coverage_summaries_destroy(summaries);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}

FC_TEST(bamcoverage_suite, test_coverage_process)
{
    const int lengths[7] = { 120, 60, 180, 70, 220, 80, 260 };
    const int32_t read_offset = 3;
    const uint32_t read_length = 10;

    char out_dir[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_cov_wXXXXXX";
    TEST_ASSERT_NOT_NULL(fc_test_mkdtemp(out_dir, sizeof(out_dir)));
    TEST_ASSERT_EQUAL_INT(0, rmdir(out_dir));

    bam_hdr_t *hdr = make_header(lengths, 7);
    cov_writer w = init_coverage_writer(
        out_dir,
        false, // per-base output off for unit tests
        true,  // write legacy coverage tables
        false, // don't inspect CIGAR op-by-op in unit tests
        -1, -1, false,
        hdr, NULL,
        NULL, NULL, 0,
        NULL, 0,
        NULL, 0, NULL);
    TEST_ASSERT_NOT_NULL(w);

    size_t prev_buf_size = 0;
    for (int32_t contig = 0; contig < 7; ++contig) {
        size_t active_len = (size_t)lengths[contig] + 1u;
        int32_t first_read_start = contig;  // do be different for each contig
        int32_t max_start = lengths[contig] - (int32_t)read_length;
        size_t n_reads = 1u + (size_t)((max_start - first_read_start) / read_offset);

        // add first read: triggers contig switch and reset of state
        {
            int32_t read_start = first_read_start;
            char qname[64];
            FC_SNPRINTF_CHECK(qname, sizeof(qname), "t%d_r0", contig);
            bam1_t *read = make_read(qname, contig, read_start, read_length, NULL, 0, NULL);
            coverage_process(w, read);
            bam_destroy1(read);
        }

        // check contig meta data
        TEST_ASSERT_EQUAL_INT(contig, w->tid);
        TEST_ASSERT_EQUAL_INT64(lengths[contig], w->contig_len);
        {
            char expected_chrom[16];
            FC_SNPRINTF_CHECK(expected_chrom, sizeof(expected_chrom), "chr%d", contig + 1);
            TEST_ASSERT_EQUAL_STRING(expected_chrom, w->chrom);
        }
        // check buffer size, contigent on if size went up or down from previous
        if (contig > 0) {
            if (active_len > prev_buf_size) {
                TEST_ASSERT_TRUE(w->buf_size > prev_buf_size);
            } else {
                TEST_ASSERT_TRUE(w->buf_size == prev_buf_size);
            }
        } else {
            TEST_ASSERT_TRUE(w->buf_size >= active_len);
        }

        // buffers are correctly zeroed on contig change:
        // after first read on this contig, only one +1/-1 pair should be present
        // we don't care about anything after active_len (memset optimisation)
        {
            int32_t *expected = calloc(active_len, sizeof(*expected));
            TEST_ASSERT_NOT_NULL(expected);
            expected[first_read_start] = 1;
            expected[first_read_start + read_length] = -1;
            TEST_ASSERT_EQUAL_INT32_ARRAY(expected, w->diff_fwd, active_len);
            TEST_ASSERT_EACH_EQUAL_INT32(0, w->diff_rev, active_len);
            free(expected);
        }

        // add remaining reads for this contig
        for (size_t read_idx = 1; read_idx < n_reads; ++read_idx) {
            int32_t read_start = first_read_start + (int32_t)read_idx * read_offset;
            char qname[64];
            FC_SNPRINTF_CHECK(qname, sizeof(qname), "t%d_r%zu", contig, read_idx);
            bam1_t *read = make_read(qname, contig, read_start, read_length, NULL, 0, NULL);
            coverage_process(w, read);
            bam_destroy1(read);
        }

        // buffers contain correct contents after all reads added for this contig
        {
            int32_t *expected_fwd = calloc(active_len, sizeof(*expected_fwd));
            TEST_ASSERT_NOT_NULL(expected_fwd);
            // trivial reimplementation of expected values
            for (size_t read_idx = 0; read_idx < n_reads; ++read_idx) {
                size_t start = (size_t)(first_read_start + (int32_t)read_idx * read_offset);
                size_t end = start + (size_t)read_length;
                TEST_ASSERT_TRUE(end < active_len);
                expected_fwd[start] += 1;
                expected_fwd[end] -= 1;
            }
            TEST_ASSERT_EQUAL_INT32_ARRAY(expected_fwd, w->diff_fwd, active_len);
            TEST_ASSERT_EACH_EQUAL_INT32(0, w->diff_rev, active_len);

            // add in a reverse read to check that diff_rev is unaffected by forward reads and vice versa
            int32_t read_start = first_read_start;
            char qname[64];
            FC_SNPRINTF_CHECK(qname, sizeof(qname), "t%d_r_rev", contig);
            uint16_t flags_rev = BAM_FREVERSE;
            bam1_t *read = make_read(qname, contig, read_start, read_length, NULL, 0, &flags_rev);
            coverage_process(w, read);
            bam_destroy1(read);
       
            TEST_ASSERT_EQUAL_INT32_ARRAY(expected_fwd, w->diff_fwd, active_len);
        
            int32_t *expected_rev = calloc(active_len, sizeof(*expected_rev));
            TEST_ASSERT_NOT_NULL(expected_rev);
            expected_rev[first_read_start] = 1;
            expected_rev[first_read_start + read_length] = -1;
            TEST_ASSERT_EQUAL_INT32_ARRAY(expected_rev, w->diff_rev, active_len);
            free(expected_rev);
            free(expected_fwd);
        }
        
        prev_buf_size = w->buf_size;
    }

    destroy_coverage_writer(w);
    sam_hdr_destroy(hdr);
    cleanup_output_dir(out_dir);
}
