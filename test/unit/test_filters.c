#include "unity.h"
#include "htslib/sam.h"

#include "fc_bam_utils.h"
#include "fc_dust.h"
#include "fc_filter.h"

static void test_filter_length(void) {
    fc_filter_opts opts = {
        .min_length = 5,
        .max_length = 10,
        .min_qscore = 0.0f
    };
    TEST_ASSERT_EQUAL_INT(FC_FILTER_TOO_SHORT, fc_filter_length(4, &opts));
    TEST_ASSERT_EQUAL_INT(FC_FILTER_OK, fc_filter_length(5, &opts));
    TEST_ASSERT_EQUAL_INT(FC_FILTER_OK, fc_filter_length(10, &opts));
    TEST_ASSERT_EQUAL_INT(FC_FILTER_TOO_LONG, fc_filter_length(11, &opts));
}

static void test_filter_quality(void) {
    fc_filter_opts opts = {
        .min_length = 0,
        .max_length = (size_t)-1,
        .min_qscore = 7.5f
    };
    TEST_ASSERT_EQUAL_INT(FC_FILTER_LOW_QUALITY, fc_filter_quality(7.4f, &opts));
    TEST_ASSERT_EQUAL_INT(FC_FILTER_OK, fc_filter_quality(7.5f, &opts));
}

static void test_filter_dust(void) {
    const uint8_t simple_seq[] =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const uint8_t complex_seq[] =
        "ACGTTGCACTGATCGTACGATCGTATGCTAGC"
        "TGACGTAGCTACGATGCTAGCATCGTACGTAG"
        "CTAGTCGATGCTAGCATGCTACGATCGTAGCT"
        "AGCATCGTACGATCGTACGATGCTAGCTAGTA";
    const int t = 20;
    const int w = 64;
    double simple_frac = fc_dust_fraction(simple_seq, sizeof(simple_seq) - 1, t, w);
    double complex_frac = fc_dust_fraction(complex_seq, sizeof(complex_seq) - 1, t, w);
    double threshold = (simple_frac + complex_frac) / 2.0;

    // simple sequence should have higher low-complexity fraction
    TEST_ASSERT_TRUE(simple_frac > complex_frac);
    // low-complexity sequence should be filtered out
    TEST_ASSERT_FALSE(fc_filter_dust(simple_seq, sizeof(simple_seq) - 1, t, w, threshold));
    // high-complexity sequence should pass
    TEST_ASSERT_TRUE(fc_filter_dust(complex_seq, sizeof(complex_seq) - 1, t, w, threshold));
}

static void test_bam_primary_flags(void) {
    bam1_t* b = bam_init1();
    TEST_ASSERT_NOT_NULL(b);

    b->core.flag = 0;
    TEST_ASSERT_TRUE(fc_is_primary(b));
    TEST_ASSERT_TRUE(fc_is_primary_qc_pass(b));

    b->core.flag = BAM_FSECONDARY;
    TEST_ASSERT_FALSE(fc_is_primary(b));

    b->core.flag = BAM_FUNMAP;
    TEST_ASSERT_FALSE(fc_is_primary(b));

    b->core.flag = BAM_FQCFAIL;
    TEST_ASSERT_TRUE(fc_is_primary(b));
    TEST_ASSERT_FALSE(fc_is_primary_qc_pass(b));

    b->core.flag = BAM_FDUP;
    TEST_ASSERT_TRUE(fc_is_primary(b));
    TEST_ASSERT_FALSE(fc_is_primary_qc_pass(b));

    bam_destroy1(b);
}

static void test_bam_tag_int_filter(void) {
    bam1_t* b = bam_init1();
    TEST_ASSERT_NOT_NULL(b);

    char tag_name[2] = {'N', 'M'};
    TEST_ASSERT_FALSE(fc_filter_int_tag(b, tag_name, 7, false));
    TEST_ASSERT_TRUE(fc_filter_int_tag(b, tag_name, 7, true));

    TEST_ASSERT_TRUE(bam_aux_update_int(b, "NM", 7) >= 0);
    TEST_ASSERT_TRUE(fc_filter_int_tag(b, tag_name, 7, false));
    TEST_ASSERT_FALSE(fc_filter_int_tag(b, tag_name, 8, false));

    bam_destroy1(b);
}

static void test_bam_rg_filter(void) {
    bam1_t* b = bam_init1();
    TEST_ASSERT_NOT_NULL(b);

    TEST_ASSERT_TRUE(fc_filter_rg(b, NULL, false));
    TEST_ASSERT_FALSE(fc_filter_rg(b, "rg1", false));
    TEST_ASSERT_TRUE(fc_filter_rg(b, "rg1", true));

    TEST_ASSERT_TRUE(bam_aux_update_str(b, "RG", 3, "rg1") >= 0);
    TEST_ASSERT_TRUE(fc_filter_rg(b, "rg1", false));
    TEST_ASSERT_FALSE(fc_filter_rg(b, "rg2", false));

    bam_destroy1(b);
}

void filter_suite(void) {
    RUN_TEST(test_filter_length);
    RUN_TEST(test_filter_quality);
    RUN_TEST(test_filter_dust);
    RUN_TEST(test_bam_primary_flags);
    RUN_TEST(test_bam_tag_int_filter);
    RUN_TEST(test_bam_rg_filter);
}
