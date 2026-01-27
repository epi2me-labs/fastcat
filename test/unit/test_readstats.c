#include <string.h>
#include <stdlib.h>

#include "unity.h"
#include "htslib/sam.h"

#include "bamstats/readstats.h"
#include "fc_read.h"

// helper to build a simple BAM record with given qualities and optional QS tag
static bam1_t* build_read(const char* qname, const char* seq, const char* qual, uint16_t flag, float qs_tag, bool add_qs) {
    int32_t l_seq = (int32_t)strlen(seq);
    int32_t l_qname = (int32_t)strlen(qname) + 1; // include null

    bam1_t* b = bam_init1();
    bam_set1(
        b, l_qname, qname,
        flag, -1, 0, 0,
        0, NULL,
        -1, 0, 0,
        l_seq, seq, qual,
        0);

    // normalise qualities to numeric phred (bam stores raw)
    uint8_t* q = bam_get_qual(b);
    for (int i = 0; i < l_seq; ++i) {
        q[i] = (uint8_t)(qual[i] - 33);
    }

    if (add_qs) {
        bam_aux_append(b, "QS", 'f', 4, (uint8_t*)&qs_tag);
    }

    return b;
}

static void test_resolve_mean_quality_prefers_qs(void) {
    bam1_t* b = build_read("r1", "ACGT", "!!!!", BAM_FUNMAP, 42.0f, true);
    read_meta meta = fc_read_meta_create(b, NULL, NULL);
    meta->mean_q = 42.0f;
    float mq = resolve_mean_quality(b, meta, false);
    fc_read_meta_destroy(meta);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 42.0f, mq);
    bam_destroy1(b);
}

static void test_resolve_mean_quality_force_recalc(void) {
    bam1_t* b = build_read("r1", "ACGT", "!!!!", BAM_FUNMAP, 42.0f, true);
    read_meta meta = fc_read_meta_create(b, NULL, NULL);
    meta->mean_q = 42.0f;
    float mq = resolve_mean_quality(b, meta, true);
    fc_read_meta_destroy(meta);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 0.0f, mq);
    bam_destroy1(b);
}

static void test_resolve_mean_quality_without_qs(void) {
    bam1_t* b = build_read("r1", "ACGT", "5555", BAM_FUNMAP, -1.0f, false);
    read_meta meta = fc_read_meta_create(b, NULL, NULL);
    meta->mean_q = -1.0f;
    float mq = resolve_mean_quality(b, meta, false);
    fc_read_meta_destroy(meta);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 20.0f, mq);
    bam_destroy1(b);
}

void readstats_suite(void) {
    RUN_TEST(test_resolve_mean_quality_prefers_qs);
    RUN_TEST(test_resolve_mean_quality_force_recalc);
    RUN_TEST(test_resolve_mean_quality_without_qs);
}
