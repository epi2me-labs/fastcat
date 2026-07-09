#include <math.h>
#include <stdio.h>
#include <limits.h>

#include "unity.h"
#include "main.h"
#include "core/fc_stats.h"
#include "core/fc_common.h"

FC_TEST(stats_suite, test_length_bins)
{
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

FC_TEST(stats_suite, test_qual_binning)
{
    read_stats* stats = create_qual_stats(0.5);
    add_qual_count(stats, 10.0);
    add_qual_count(stats, 10.4);
    add_qual_count(stats, 100.0);   // on boundary
    add_qual_count(stats, 120.0);   // clamped

    TEST_ASSERT_EQUAL_size_t(2, stats->counts[20]); // 10.0 and 10.4 in bucket 20 (10.0 / 0.5)
    TEST_ASSERT_EQUAL_size_t(2, stats->counts[200]); // 100 and clamped 120 go to bucket 200

    destroy_qual_stats(stats);
}

FC_TEST(stats_suite, test_leading_decimals)
{
    TEST_ASSERT_EQUAL_size_t(7, _leading_decimals(1.0));
    TEST_ASSERT_EQUAL_size_t(1, _leading_decimals(0.1));
    TEST_ASSERT_EQUAL_size_t(2, _leading_decimals(0.01));
    TEST_ASSERT_EQUAL_size_t(3, _leading_decimals(0.001));
}

FC_TEST(stats_suite, test_mean_qual_naive)
{
    char q[] = {'5', '5', '5', '5'}; // phred 20
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 20.0f, mean_qual_naive(q, 4));
}

FC_TEST(stats_suite, test_length_stats_layout)
{
    read_stats* stats = create_length_stats();
    TEST_ASSERT_EQUAL_size_t(10000001, stats->n);
    TEST_ASSERT_EQUAL_size_t(0, stats->edges[0]);
    TEST_ASSERT_EQUAL_size_t(1, stats->edges[1]);
    TEST_ASSERT_EQUAL_size_t(10000000, stats->edges[stats->n - 1]);
    TEST_ASSERT_EQUAL_size_t(SIZE_MAX, stats->edges[stats->n]);
    // spot check monotonicity
    TEST_ASSERT_TRUE(stats->edges[5000] < stats->edges[5001]);
    destroy_length_stats(stats);
}

FC_TEST(stats_suite, test_print_stats_length_overflow_bucket_output)
{
    read_stats *stats = create_length_stats();
    add_length_count(stats, 10000001);

    FILE *fp = tmpfile();
    TEST_ASSERT_NOT_NULL(fp);
    print_stats(stats, false, true, fp);
    TEST_ASSERT_EQUAL_INT(0, fflush(fp));
    rewind(fp);

    char line[128] = {0};
    TEST_ASSERT_NOT_NULL(fgets(line, sizeof(line), fp));
    char expected[128] = {0};
    FC_SNPRINTF_CHECK(expected, sizeof(expected), "10000000\t%zu\t1\n", (size_t)SIZE_MAX);
    TEST_ASSERT_EQUAL_STRING(expected, line);
    TEST_ASSERT_NULL(fgets(line, sizeof(line), fp));

    fclose(fp);
    destroy_length_stats(stats);
}

FC_TEST(stats_suite, test_print_stats_tsv_output)
{
    read_stats *stats = create_qual_stats(0.5f);
    add_qual_count(stats, 10.0f);
    add_qual_count(stats, 10.4f);

    FILE *fp = tmpfile();
    TEST_ASSERT_NOT_NULL(fp);
    print_stats(stats, false, true, fp);
    TEST_ASSERT_EQUAL_INT(0, fflush(fp));
    rewind(fp);

    char line[128] = {0};
    TEST_ASSERT_NOT_NULL(fgets(line, sizeof(line), fp));
    TEST_ASSERT_EQUAL_STRING("10.0\t10.5\t2\n", line);
    TEST_ASSERT_NULL(fgets(line, sizeof(line), fp));

    fclose(fp);
    destroy_qual_stats(stats);
}

FC_TEST(stats_suite, test_print_stats_pretty_output)
{
    read_stats *stats = create_qual_stats(0.5f);
    add_qual_count(stats, 10.0f);
    add_qual_count(stats, 10.4f);

    FILE *fp = tmpfile();
    TEST_ASSERT_NOT_NULL(fp);
    print_stats(stats, false, false, fp);
    TEST_ASSERT_EQUAL_INT(0, fflush(fp));
    rewind(fp);

    char line[128] = {0};
    TEST_ASSERT_NOT_NULL(fgets(line, sizeof(line), fp));
    TEST_ASSERT_EQUAL_STRING("[10.0, 10.5)\t2\n", line);
    TEST_ASSERT_NULL(fgets(line, sizeof(line), fp));

    fclose(fp);
    destroy_qual_stats(stats);
}
