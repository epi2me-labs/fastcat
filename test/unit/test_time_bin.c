#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "htslib/sam.h"

#include "unity.h"
#include "core/fc_time_bin.h"

static void assert_parse_time(const char *text, double expected)
{
    double observed = NAN;
    TEST_ASSERT_TRUE(fc_time_bin_parse_timestamp(text, &observed));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, expected, observed);
}

static double row_total_reads(const fc_time_bin_table *table)
{
    double total = 0.0;
    for (size_t i = 0; i < table->n; i++) {
        total += (double)table->rows[i].read_count;
    }
    return total;
}

static fc_read *make_time_bin_test_read(const char *name)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);

    const char *seq = "ACGT";
    const char *qual_ascii = "IIII";
    TEST_ASSERT_TRUE(bam_set1(
        read->bam,
        strlen(name) + 1,
        name,
        BAM_FUNMAP,
        -1,
        -1,
        0,
        0,
        NULL,
        -1,
        -1,
        0,
        strlen(seq),
        seq,
        qual_ascii,
        0) >= 0);
    return read;
}

static void add_time_bin_test_read(
    fc_time_bin *model,
    size_t barcode,
    const char *timestamp,
    size_t length,
    double mean_q)
{
    bam1_t bam = { 0 };
    read_meta meta = {
        .start_time = (char *)timestamp,
        .mean_q = mean_q
    };
    fc_read read = {
        .bam = &bam,
        .meta = &meta
    };
    bam.core.l_qseq = (int32_t)length;
    fc_time_bin_add_read(model, barcode, &read);
}

static void test_time_bin_parse_timestamp_variants(void)
{
    assert_parse_time("123.5", 123.5);
    assert_parse_time("1970-01-01T00:00:00Z", 0.0);
    assert_parse_time("1970-01-01T00:00:00.5Z", 0.0);
    assert_parse_time("1970-01-01T01:00:00+01:00", 3600.0);
    assert_parse_time("1970-01-01T00:30:00-0030", 1800.0);
    assert_parse_time("2024-01-24T17:15:41.615+00:00", 1706116541.0);
    assert_parse_time("2025-01-13T10:45:28.681306+00:00", 1736765128.0);
    assert_parse_time("2016-01-19T15:21:32.59+02:00", 1453216892.0);

    double out = 0.0;
    TEST_ASSERT_FALSE(fc_time_bin_parse_timestamp("", &out));
    TEST_ASSERT_FALSE(fc_time_bin_parse_timestamp("not-a-time", &out));
}

static void test_time_bin_counts_means_quantiles_and_cumulative_totals(void)
{
    fc_time_bin *model = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    add_time_bin_test_read(model, 0, "100", 100, 10.0);
    add_time_bin_test_read(model, 0, "105", 200, 20.0);
    add_time_bin_test_read(model, 0, "110", 300, 30.0);

    fc_time_bin_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(2, table.n);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, row_total_reads(&table));

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, table.rows[0].time_bin_start_s);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 10.0, table.rows[0].time_bin_end_s);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(300, table.rows[0].base_count);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(300, table.rows[0].cumulative_base_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 150.0, table.rows[0].length_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 105.0, table.rows[0].length_quantiles.q05);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 125.0, table.rows[0].length_quantiles.q25);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 150.0, table.rows[0].length_quantiles.q50);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 175.0, table.rows[0].length_quantiles.q75);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 195.0, table.rows[0].length_quantiles.q95);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 15.0, table.rows[0].qscore_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 10.5, table.rows[0].qscore_quantiles.q05);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 12.5, table.rows[0].qscore_quantiles.q25);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 15.0, table.rows[0].qscore_quantiles.q50);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 17.5, table.rows[0].qscore_quantiles.q75);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 19.5, table.rows[0].qscore_quantiles.q95);

    TEST_ASSERT_EQUAL_size_t(1, table.rows[1].read_count);
    TEST_ASSERT_EQUAL_size_t(300, table.rows[1].base_count);
    TEST_ASSERT_EQUAL_size_t(3, table.rows[1].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(600, table.rows[1].cumulative_base_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 300.0, table.rows[1].length_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 30.0, table.rows[1].qscore_mean);

    fc_time_bin_table_destroy(&table);
    fc_time_bin_destroy(model);
}

static void test_time_bin_dynamic_backward_growth_preserves_empty_bins(void)
{
    fc_time_bin *model = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    add_time_bin_test_read(model, 0, "100", 100, 10.0);
    add_time_bin_test_read(model, 0, "80", 80, 8.0);

    fc_time_bin_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(3, table.n);

    TEST_ASSERT_EQUAL_size_t(1, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(80, table.rows[0].base_count);
    TEST_ASSERT_EQUAL_size_t(1, table.rows[0].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(80, table.rows[0].cumulative_base_count);

    TEST_ASSERT_EQUAL_size_t(0, table.rows[1].read_count);
    TEST_ASSERT_TRUE(isnan(table.rows[1].length_mean));
    TEST_ASSERT_TRUE(isnan(table.rows[1].length_quantiles.q05));
    TEST_ASSERT_TRUE(isnan(table.rows[1].length_quantiles.q95));
    TEST_ASSERT_TRUE(isnan(table.rows[1].qscore_mean));
    TEST_ASSERT_TRUE(isnan(table.rows[1].qscore_quantiles.q05));
    TEST_ASSERT_TRUE(isnan(table.rows[1].qscore_quantiles.q95));
    TEST_ASSERT_EQUAL_size_t(1, table.rows[1].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(80, table.rows[1].cumulative_base_count);

    TEST_ASSERT_EQUAL_size_t(1, table.rows[2].read_count);
    TEST_ASSERT_EQUAL_size_t(100, table.rows[2].base_count);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[2].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(180, table.rows[2].cumulative_base_count);

    fc_time_bin_table_destroy(&table);
    fc_time_bin_destroy(model);
}

static void test_time_bin_materialize_merged_limits_rows_and_merges_sketches(void)
{
    fc_time_bin *model = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    char timestamp[32];
    for (size_t i = 0; i < 48; i++) {
        snprintf(timestamp, sizeof(timestamp), "%zu", i * 10);
        add_time_bin_test_read(model, 0, timestamp, 100 + i, 10.0 + (double)i);
    }

    fc_time_bin_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize_merged(model, 0, 24, &table));
    TEST_ASSERT_EQUAL_size_t(24, table.n);

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, table.rows[0].time_bin_start_s);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 20.0, table.rows[0].time_bin_end_s);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(201, table.rows[0].base_count);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(201, table.rows[0].cumulative_base_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 100.5, table.rows[0].length_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 100.5, table.rows[0].length_quantiles.q50);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 10.5, table.rows[0].qscore_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 10.5, table.rows[0].qscore_quantiles.q50);

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 460.0, table.rows[23].time_bin_start_s);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 480.0, table.rows[23].time_bin_end_s);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[23].read_count);
    TEST_ASSERT_EQUAL_size_t(293, table.rows[23].base_count);
    TEST_ASSERT_EQUAL_size_t(48, table.rows[23].cumulative_read_count);
    TEST_ASSERT_EQUAL_size_t(5928, table.rows[23].cumulative_base_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 146.5, table.rows[23].length_quantiles.q50);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 56.5, table.rows[23].qscore_quantiles.q50);

    fc_time_bin_table_destroy(&table);
    fc_time_bin_destroy(model);
}

static void test_time_bin_barcode_state_is_isolated(void)
{
    fc_time_bin *model = fc_time_bin_create(3, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    add_time_bin_test_read(model, 1, "0", 100, 10.0);
    add_time_bin_test_read(model, 2, "0", 200, 20.0);

    TEST_ASSERT_FALSE(fc_time_bin_seen(model, 0));
    TEST_ASSERT_TRUE(fc_time_bin_seen(model, 1));
    TEST_ASSERT_TRUE(fc_time_bin_seen(model, 2));
    TEST_ASSERT_FALSE(fc_time_bin_seen(model, 3));

    fc_time_bin_table table0 = { 0 };
    fc_time_bin_table table1 = { 0 };
    fc_time_bin_table table2 = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 0, &table0));
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 1, &table1));
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 2, &table2));
    TEST_ASSERT_EQUAL_INT(-1, fc_time_bin_materialize(model, 3, &table2));

    TEST_ASSERT_EQUAL_size_t(0, table0.n);
    TEST_ASSERT_EQUAL_size_t(1, table1.n);
    TEST_ASSERT_EQUAL_size_t(1, table2.n);
    TEST_ASSERT_EQUAL_size_t(100, table1.rows[0].base_count);
    TEST_ASSERT_EQUAL_size_t(200, table2.rows[0].base_count);

    fc_time_bin_table_destroy(&table0);
    fc_time_bin_table_destroy(&table1);
    fc_time_bin_table_destroy(&table2);
    fc_time_bin_destroy(model);
}

static void test_time_bin_missing_and_invalid_timestamps_are_counted(void)
{
    fc_time_bin *model = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    read_meta meta = { 0 };
    fc_read read = {
        .bam = bam,
        .meta = &meta
    };

    meta.start_time = "";
    fc_time_bin_add_read(model, 0, &read);
    meta.start_time = "not-a-time";
    fc_time_bin_add_read(model, 0, &read);

    TEST_ASSERT_TRUE(fc_time_bin_seen(model, 0));
    fc_time_bin_skip_counts skipped = fc_time_bin_skipped_counts(model);
    TEST_ASSERT_EQUAL_size_t(1, skipped.missing);
    TEST_ASSERT_EQUAL_size_t(1, skipped.invalid);

    fc_time_bin_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(0, table.n);

    fc_time_bin_table_destroy(&table);
    bam_destroy1(bam);
    fc_time_bin_destroy(model);
}

static void test_time_bin_add_read_uses_bam_st_z_timestamp_tags(void)
{
    fc_time_bin *model = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(model);

    fc_read *read1 = make_time_bin_test_read("r1");
    const char *time1 = "1970-01-01T00:01:40Z";
    TEST_ASSERT_TRUE(bam_aux_update_str(read1->bam, "st", (int)strlen(time1), time1) >= 0);
    TEST_ASSERT_TRUE(bam_aux_update_float(read1->bam, "qs", 10.0f) >= 0);
    fc_read_reset_from_bam(read1, NULL, NULL);
    TEST_ASSERT_EQUAL_STRING(time1, read1->meta->start_time);
    fc_time_bin_add_read(model, 0, read1);

    fc_read *read2 = make_time_bin_test_read("r2");
    const char *time2 = "1970-01-01T00:01:45.5Z";
    TEST_ASSERT_TRUE(bam_aux_update_str(read2->bam, "st", (int)strlen(time2), time2) >= 0);
    TEST_ASSERT_TRUE(bam_aux_update_float(read2->bam, "qs", 20.0f) >= 0);
    fc_read_reset_from_bam(read2, NULL, NULL);
    double parsed = 0.0;
    TEST_ASSERT_TRUE(fc_time_bin_parse_timestamp(read2->meta->start_time, &parsed));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 105.0, parsed);
    fc_time_bin_add_read(model, 0, read2);

    fc_time_bin_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_time_bin_materialize(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(1, table.n);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(8, table.rows[0].base_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 4.0, table.rows[0].length_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 15.0, table.rows[0].qscore_mean);

    fc_time_bin_table_destroy(&table);
    fc_read_destroy(read1);
    fc_read_destroy(read2);
    fc_time_bin_destroy(model);
}

void time_bin_suite(void)
{
    RUN_TEST(test_time_bin_parse_timestamp_variants);
    RUN_TEST(test_time_bin_counts_means_quantiles_and_cumulative_totals);
    RUN_TEST(test_time_bin_dynamic_backward_growth_preserves_empty_bins);
    RUN_TEST(test_time_bin_materialize_merged_limits_rows_and_merges_sketches);
    RUN_TEST(test_time_bin_barcode_state_is_isolated);
    RUN_TEST(test_time_bin_missing_and_invalid_timestamps_are_counted);
    RUN_TEST(test_time_bin_add_read_uses_bam_st_z_timestamp_tags);
}
