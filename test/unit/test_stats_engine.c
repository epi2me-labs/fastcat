#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "htslib/sam.h"

#include "fc_stats_engine.h"

static size_t sum_counts(const read_stats* stats) {
    size_t total = 0;
    if (stats == NULL) return 0;
    for (size_t i = 0; i < stats->n; ++i) {
        total += stats->counts[i];
    }
    return total;
}

static bam1_t* build_unmapped_read(const char* qname, const char* seq, const char* qual) {
    int32_t l_seq = (int32_t)strlen(seq);
    int32_t l_qname = (int32_t)strlen(qname) + 1;

    bam1_t* b = bam_init1();
    bam_set1(
        b, l_qname, qname,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        l_seq, seq, qual,
        0);

    uint8_t* q = bam_get_qual(b);
    for (int i = 0; i < l_seq; ++i) {
        q[i] = (uint8_t)(qual[i] - 33);
    }

    return b;
}

static bam1_t* build_mapped_read(const char* qname, const char* seq, const char* qual, int32_t tid, int32_t pos) {
    int32_t l_seq = (int32_t)strlen(seq);
    int32_t l_qname = (int32_t)strlen(qname) + 1;
    uint32_t cigar = bam_cigar_gen(l_seq, BAM_CMATCH);

    bam1_t* b = bam_init1();
    bam_set1(
        b, l_qname, qname,
        0, tid, pos, 60,
        1, &cigar,
        -1, -1, 0,
        l_seq, seq, qual,
        0);

    uint8_t* q = bam_get_qual(b);
    for (int i = 0; i < l_seq; ++i) {
        q[i] = (uint8_t)(qual[i] - 33);
    }
    TEST_ASSERT_TRUE(bam_aux_update_int(b, "NM", 0) >= 0);

    return b;
}

static void test_stats_engine_unmapped_tracked_in_unmapped_stats(void) {
    bam1_t* b = build_unmapped_read("r1", "ACGT", "!!!!");
    fc_read* fr = fc_read_from_bam(b, NULL, false);

    fc_stats_engine* engine = fc_stats_engine_create(true, false, false, 0.0f, 0.0f, false);
    TEST_ASSERT_EQUAL_INT(FC_STATS_UNMAPPED, fc_stats_engine_add_read(engine, fr, NULL, NULL));
    TEST_ASSERT_EQUAL_UINT64(0, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->length_stats_unmapped));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(fr);
    bam_destroy1(b);
}

static void test_stats_engine_unmapped_ignored_without_tracking(void) {
    bam1_t* b = build_unmapped_read("r1", "ACGT", "!!!!");
    fc_read* fr = fc_read_from_bam(b, NULL, false);

    fc_stats_engine* engine = fc_stats_engine_create(false, false, false, 0.0f, 0.0f, false);
    TEST_ASSERT_EQUAL_INT(FC_STATS_UNMAPPED, fc_stats_engine_add_read(engine, fr, NULL, NULL));
    TEST_ASSERT_EQUAL_UINT64(0, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(0, sum_counts(engine->length_stats_unmapped));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(fr);
    bam_destroy1(b);
}

static void test_stats_engine_aligned_metrics(void) {
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "100", NULL);

    bam1_t* b = build_mapped_read("r1", "ACGTACGTAC", "!!!!!!!!!!", 0, 5);
    fc_read* fr = fc_read_from_bam(b, hdr, false);

    fc_stats_engine* engine = fc_stats_engine_create(false, true, false, 0.0f, 0.0f, false);
    fc_stats_aligned metrics = {0};
    TEST_ASSERT_EQUAL_INT(FC_STATS_ALIGNED, fc_stats_engine_add_read(engine, fr, hdr, &metrics));

    TEST_ASSERT_EQUAL_UINT32(10, metrics.read_length);
    TEST_ASSERT_EQUAL_UINT32(5, metrics.rstart);
    TEST_ASSERT_EQUAL_UINT32(15, metrics.rend);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.aligned_ref_len);
    TEST_ASSERT_EQUAL_UINT32(100, metrics.ref_length);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.qstart);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.qend);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.match);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.ins);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.delt);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.sub);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.length);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.coverage);
    TEST_ASSERT_EQUAL_FLOAT(10.0f, metrics.ref_cover);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.iden);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.acc);
    TEST_ASSERT_EQUAL_CHAR('+', metrics.direction);

    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->qual_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->acc_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->cov_stats));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(fr);
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
}

void stats_engine_suite(void) {
    RUN_TEST(test_stats_engine_unmapped_tracked_in_unmapped_stats);
    RUN_TEST(test_stats_engine_unmapped_ignored_without_tracking);
    RUN_TEST(test_stats_engine_aligned_metrics);
}
