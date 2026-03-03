#include <math.h>
#include <stdint.h>

#include "unity.h"
#include "core/fc_common.h"

static void test_mean_qual_from_bam_matches_naive(void) {
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.0f, mean_qual_from_bam_naive((uint8_t[]){20, 20, 20, 20}, 4));
    TEST_ASSERT_FLOAT_WITHIN(1e-6,
        mean_qual_from_bam_naive((uint8_t[]){20, 20, 20, 20}, 4),
        mean_qual_from_bam((uint8_t[]){20, 20, 20, 20}, 4));
}

static void test_mean_qual_from_bam_sentinel_nan(void) {
    uint8_t q[] = {0xff, 10, 10};
    float v = mean_qual_from_bam(q, 3);
    TEST_ASSERT_TRUE(isnan(v));
}

static void test_mean_qual_empty_ascii(void) {
    char* q = NULL;
    float v = mean_qual(q, 0);
    TEST_ASSERT_TRUE(isnan(v));
}

static void test_mean_qual_empty_bam(void) {
    uint8_t* q = NULL;
    float v = mean_qual_from_bam(q, 0);
    TEST_ASSERT_TRUE(isnan(v));
}

static void test_mean_qual_from_bam_high_phred_bound(void) {
    uint8_t q[] = {99, 99, 99, 99};
    float fast = mean_qual_from_bam(q, 4);
    float naive = mean_qual_from_bam_naive(q, 4);
    TEST_ASSERT_FLOAT_WITHIN(1e-5, naive, fast);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 99.0f, fast);
}

void quality_suite(void) {
    RUN_TEST(test_mean_qual_from_bam_matches_naive);
    RUN_TEST(test_mean_qual_from_bam_sentinel_nan);
    RUN_TEST(test_mean_qual_empty_ascii);
    RUN_TEST(test_mean_qual_empty_bam);
    RUN_TEST(test_mean_qual_from_bam_high_phred_bound);
}
