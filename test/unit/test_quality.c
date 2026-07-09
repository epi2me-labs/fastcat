#include <math.h>
#include <stdint.h>

#include "unity.h"
#include "main.h"
#include "core/fc_common.h"

FC_TEST(quality_suite, test_mean_qual_from_bam_matches_naive)
{
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.0f, mean_qual_from_bam_naive((uint8_t[]){20, 20, 20, 20}, 4));
    TEST_ASSERT_FLOAT_WITHIN(1e-6,
        mean_qual_from_bam_naive((uint8_t[]){20, 20, 20, 20}, 4),
        mean_qual_from_bam((uint8_t[]){20, 20, 20, 20}, 4));
}

FC_TEST(quality_suite, test_mean_qual_from_bam_sentinel_nan)
{
    uint8_t q[] = {0xff, 10, 10};
    float v = mean_qual_from_bam(q, 3);
    TEST_ASSERT_TRUE(isnan(v));
}

FC_TEST(quality_suite, test_mean_qual_empty_ascii)
{
    char* q = NULL;
    float v = mean_qual(q, 0);
    TEST_ASSERT_TRUE(isnan(v));
}

FC_TEST(quality_suite, test_mean_qual_empty_bam)
{
    uint8_t* q = NULL;
    float v = mean_qual_from_bam(q, 0);
    TEST_ASSERT_TRUE(isnan(v));
}

FC_TEST(quality_suite, test_mean_qual_from_bam_high_phred_bound)
{
    uint8_t q[] = {99, 99, 99, 99};
    float fast = mean_qual_from_bam(q, 4);
    float naive = mean_qual_from_bam_naive(q, 4);
    TEST_ASSERT_FLOAT_WITHIN(1e-5, naive, fast);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 99.0f, fast);
}
