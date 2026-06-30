#include <math.h>
#include <stdio.h>

#include "unity.h"
#include "core/fc_quantile.h"
#include "core/fc_stats.h"

static void test_empty_histogram(void)
{
    read_stats *stats = create_qual_stats(1.0f);

    double q50 = fc_compute_quantile(stats, 0.5f);
    TEST_ASSERT_TRUE(isnan(q50));

    fc_quantile_triple triple = fc_compute_quantile_triple(stats);
    TEST_ASSERT_TRUE(isnan(triple.q25));
    TEST_ASSERT_TRUE(isnan(triple.q50));
    TEST_ASSERT_TRUE(isnan(triple.q75));

    destroy_qual_stats(stats);
}

static void test_uniform_distribution(void)
{
    // Create histogram with 1-unit bins from 0 to 100
    read_stats *stats = create_qual_stats(1.0f);

    // Add uniform distribution: 10 counts in each bin from 0 to 99
    for (size_t i = 0; i < 100; i++) {
        for (int j = 0; j < 10; j++) {
            add_qual_count(stats, (float)i + 0.5f); // mid-bin values
        }
    }

    // Total: 1000 values uniformly distributed
    // Q25 should be around bin 25, Q50 around 50, Q75 around 75
    fc_quantile_triple triple = fc_compute_quantile_triple(stats);

    TEST_ASSERT_FLOAT_WITHIN(1.0f, 25.0f, triple.q25);
    TEST_ASSERT_FLOAT_WITHIN(1.0f, 50.0f, triple.q50);
    TEST_ASSERT_FLOAT_WITHIN(1.0f, 75.0f, triple.q75);

    destroy_qual_stats(stats);
}

static void test_single_bin(void)
{
    read_stats *stats = create_qual_stats(1.0f);

    // Add 100 values all in bin 42: [42, 43)
    for (int i = 0; i < 100; i++) {
        add_qual_count(stats, 42.5f);
    }

    fc_quantile_triple triple = fc_compute_quantile_triple(stats);

    // With uniform distribution assumption within bin:
    // Q25 at position 25/100 = 0.25 through bin -> 42 + 0.25 = 42.25
    // Q50 at position 50/100 = 0.50 through bin -> 42 + 0.50 = 42.50
    // Q75 at position 75/100 = 0.75 through bin -> 42 + 0.75 = 42.75
    TEST_ASSERT_FLOAT_WITHIN(0.01f, 42.25f, triple.q25);
    TEST_ASSERT_FLOAT_WITHIN(0.01f, 42.50f, triple.q50);
    TEST_ASSERT_FLOAT_WITHIN(0.01f, 42.75f, triple.q75);

    destroy_qual_stats(stats);
}

static void test_zero_quantile_with_leading_empty_bins(void)
{
    read_stats *stats = create_qual_stats(1.0f);

    // Leave bins [0, 42) empty and put all observations in [42, 43).
    for (int i = 0; i < 100; i++) {
        add_qual_count(stats, 42.5f);
    }

    double q0 = fc_compute_quantile(stats, 0.0f);

    // The minimum observed bucket starts at 42, so p=0 should resolve there.
    TEST_ASSERT_FALSE(isnan(q0));
    TEST_ASSERT_FLOAT_WITHIN(0.01f, 42.0f, q0);

    destroy_qual_stats(stats);
}

static void test_skewed_distribution(void)
{
    read_stats *stats = create_qual_stats(1.0f);

    // Exponential decay: many low values, few high values
    for (int i = 0; i < 100; i++) add_qual_count(stats, 0.5f);
    for (int i = 0; i < 50; i++) add_qual_count(stats, 1.5f);
    for (int i = 0; i < 25; i++) add_qual_count(stats, 2.5f);
    for (int i = 0; i < 13; i++) add_qual_count(stats, 3.5f);
    for (int i = 0; i < 6; i++) add_qual_count(stats, 4.5f);
    for (int i = 0; i < 3; i++) add_qual_count(stats, 5.5f);
    for (int i = 0; i < 2; i++) add_qual_count(stats, 6.5f);

    fc_quantile_triple triple = fc_compute_quantile_triple(stats);

    // Verify values are reasonable and ordered
    TEST_ASSERT_FALSE(isnan(triple.q25));
    TEST_ASSERT_FALSE(isnan(triple.q50));
    TEST_ASSERT_FALSE(isnan(triple.q75));
    TEST_ASSERT_TRUE(triple.q25 < triple.q50);
    TEST_ASSERT_TRUE(triple.q50 < triple.q75);
    TEST_ASSERT_TRUE(triple.q25 < 5.0f); // heavily skewed

    destroy_qual_stats(stats);
}

static void test_segmented_bins(void)
{
    read_stats *stats = create_qc_length_bins();

    // Add some reads at different length scales
    // Short reads around 1kb
    for (int i = 0; i < 100; i++) {
        add_length_count(stats, 1000 + i);
    }

    // Medium reads around 10kb
    for (int i = 0; i < 50; i++) {
        add_length_count(stats, 10000 + i * 10);
    }

    // Long reads around 50kb
    for (int i = 0; i < 20; i++) {
        add_length_count(stats, 50000 + i * 100);
    }

    fc_quantile_triple triple = fc_compute_quantile_triple(stats);

    // Verify quantiles are reasonable
    TEST_ASSERT_FALSE(isnan(triple.q25));
    TEST_ASSERT_FALSE(isnan(triple.q50));
    TEST_ASSERT_FALSE(isnan(triple.q75));
    TEST_ASSERT_TRUE(triple.q25 < triple.q50);
    TEST_ASSERT_TRUE(triple.q50 < triple.q75);

    // Q25 should be in lower group (around 1kb)
    // Q75 should be higher
    TEST_ASSERT_TRUE(triple.q25 < 5000.0f);
    TEST_ASSERT_TRUE(triple.q75 > 10000.0f);

    destroy_length_stats(stats);
}

static void test_overflow_bin(void)
{
    // Test that overflow bin (unbounded upper edge) returns lower bound
    read_stats *stats = create_qc_length_bins();

    // Add 70 reads in normal bins (below 100kb)
    for (int i = 0; i < 70; i++) {
        add_length_count(stats, 1000 + i * 100);
    }

    // Add 30 reads in overflow bin (>= 100kb)
    for (int i = 0; i < 30; i++) {
        add_length_count(stats, 100000 + i * 1000); // 100kb to 129kb
    }

    // Total: 100 reads
    // First 70 reads: 1000-7900
    // Last 30 reads: 100000-129000
    // Q75 at position 75 falls in the overflow bin

    fc_quantile_triple triple = fc_compute_quantile_triple(stats);

    // Q25 and Q50 should be in normal bins
    TEST_ASSERT_TRUE(triple.q25 < 100000.0f);
    TEST_ASSERT_TRUE(triple.q50 < 100000.0f);

    // Q75 should fall in overflow bin (unbounded upper edge)
    // Should return 100000 (lower bound), not try to interpolate to infinity
    TEST_ASSERT_FLOAT_WITHIN(1.0f, 100000.0f, triple.q75);

    destroy_length_stats(stats);
}

static void test_fixed_width_terminal_bucket_returns_lower_bound(void)
{
    read_stats *stats = create_qual_stats(1.0f);

    // Values >= 100 are all clamped into the terminal bucket.
    for (int i = 0; i < 100; i++) {
        add_qual_count(stats, 100.0f);
    }

    double q50 = fc_compute_quantile(stats, 0.5f);

    // Once the terminal bucket is acting as overflow, quantile readback should
    // not invent values beyond the capped histogram range.
    TEST_ASSERT_FLOAT_WITHIN(0.01f, 100.0f, q50);

    destroy_qual_stats(stats);
}

static void test_fixed_width_vs_segmented(void)
{
    // Fixed-width
    read_stats *fixed = create_qual_stats(0.5f);
    for (int i = 0; i < 100; i++) {
        add_qual_count(fixed, 10.0f + (float)i * 0.1f);
    }
    fc_quantile_triple triple_fixed = fc_compute_quantile_triple(fixed);
    TEST_ASSERT_FALSE(isnan(triple_fixed.q50));

    // Segmented
    read_stats *segmented = create_qc_length_bins();
    for (int i = 0; i < 100; i++) {
        add_length_count(segmented, 1000 + i * 10);
    }
    fc_quantile_triple triple_seg = fc_compute_quantile_triple(segmented);
    TEST_ASSERT_FALSE(isnan(triple_seg.q50));

    destroy_qual_stats(fixed);
    destroy_length_stats(segmented);
}

void quantile_suite(void)
{
    RUN_TEST(test_empty_histogram);
    RUN_TEST(test_uniform_distribution);
    RUN_TEST(test_single_bin);
    RUN_TEST(test_zero_quantile_with_leading_empty_bins);
    RUN_TEST(test_skewed_distribution);
    RUN_TEST(test_segmented_bins);
    RUN_TEST(test_overflow_bin);
    RUN_TEST(test_fixed_width_terminal_bucket_returns_lower_bound);
    RUN_TEST(test_fixed_width_vs_segmented);
}
