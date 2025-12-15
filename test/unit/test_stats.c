#include <math.h>
#include <stdio.h>

#include "unity.h"
#include "stats.h"
#include "common.h"

static void test_length_bins(void) {
    read_stats* stats = create_length_stats();
    add_length_count(stats, 0);         // first bin
    add_length_count(stats, 5);         // distinct bin
    add_length_count(stats, 9999999);   // last in fixed range
    add_length_count(stats, 15000000);  // overflow bucket

    TEST_ASSERT_EQUAL_size_t(1, stats->counts[0]);
    TEST_ASSERT_EQUAL_size_t(1, stats->counts[5]);
    TEST_ASSERT_EQUAL_size_t(1, stats->counts[9999999]);
    TEST_ASSERT_EQUAL_size_t(1, stats->counts[stats->n - 1]); // overflow

    destroy_length_stats(stats);
}

static void test_qual_binning(void) {
    read_stats* stats = create_qual_stats(0.5);
    add_qual_count(stats, 10.0);
    add_qual_count(stats, 10.4);
    add_qual_count(stats, 100.0);   // on boundary
    add_qual_count(stats, 120.0);   // clamped

    TEST_ASSERT_EQUAL_size_t(2, stats->counts[20]); // 10.0 and 10.4 in bucket 20 (10.0 / 0.5)
    TEST_ASSERT_EQUAL_size_t(2, stats->counts[200]); // 100 and clamped 120 go to bucket 200

    destroy_qual_stats(stats);
}

static void test_leading_decimals(void) {
    TEST_ASSERT_EQUAL_size_t(7, _leading_decimals(1.0));
    TEST_ASSERT_EQUAL_size_t(1, _leading_decimals(0.1));
    TEST_ASSERT_EQUAL_size_t(2, _leading_decimals(0.01));
    TEST_ASSERT_EQUAL_size_t(3, _leading_decimals(0.001));
}

static void test_mean_qual_naive(void) {
    char q[] = {'5', '5', '5', '5'}; // phred 20
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 20.0f, mean_qual_naive(q, 4));
}

static void test_length_stats_layout(void) {
    read_stats* stats = create_length_stats();
    TEST_ASSERT_EQUAL_size_t(10000001, stats->n);
    TEST_ASSERT_EQUAL_size_t(0, stats->edges[0]);
    TEST_ASSERT_EQUAL_size_t(1, stats->edges[1]);
    TEST_ASSERT_EQUAL_size_t(10000000, stats->edges[stats->n - 1]);
    // spot check monotonicity
    TEST_ASSERT_TRUE(stats->edges[5000] < stats->edges[5001]);
    destroy_length_stats(stats);
}

void stats_suite(void) {
    RUN_TEST(test_length_bins);
    RUN_TEST(test_qual_binning);
    RUN_TEST(test_leading_decimals);
    RUN_TEST(test_mean_qual_naive);
    RUN_TEST(test_length_stats_layout);
}
