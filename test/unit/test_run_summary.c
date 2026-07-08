#include <math.h>
#include <stdint.h>

#include "unity.h"

#include "core/fc_run_summary.h"

static void add_run_summary_read(
    fc_run_summary *model,
    size_t barcode,
    size_t length,
    const char *sample_id,
    const char *run_id,
    const char *basecaller_model,
    const char *flow_cell_id,
    const char *barcode_alias,
    const char *start_time)
{
    bam1_t bam = { 0 };
    read_meta meta = {
        .sample_id = (char *)(sample_id ? sample_id : ""),
        .runid = (char *)(run_id ? run_id : ""),
        .basecaller = (char *)(basecaller_model ? basecaller_model : ""),
        .flow_cell_id = (char *)(flow_cell_id ? flow_cell_id : ""),
        .barcode_alias = (char *)(barcode_alias ? barcode_alias : ""),
        .start_time = (char *)(start_time ? start_time : "")
    };
    fc_read read = {
        .bam = &bam,
        .meta = &meta
    };

    bam.core.l_qseq = (int32_t)length;
    fc_run_summary_add_read(model, barcode, &read);
}

static void test_run_summary_collapses_multiple_values(void)
{
    fc_run_summary *model = fc_run_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_run_summary_read(
        model, 0, 100,
        "sample-a",
        "run-a",
        "model-a",
        "flow-a",
        "alias-a",
        "2026-01-01T00:00:00Z");
    add_run_summary_read(
        model, 0, 200,
        "sample-b",
        "run-b",
        "model-b",
        "flow-b",
        "alias-b",
        "2026-01-01T01:00:00Z");

    TEST_ASSERT_TRUE(fc_run_summary_seen(model, 0));

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, NULL, "fallback-sample", "barcode0001", &row));

    TEST_ASSERT_EQUAL_STRING("(multiple)", row.sample_id);
    TEST_ASSERT_EQUAL_STRING("barcode0001", row.barcode_id);
    TEST_ASSERT_EQUAL_STRING("alias-a", row.barcode_alias);
    TEST_ASSERT_EQUAL_STRING("(multiple)", row.run_id);
    TEST_ASSERT_EQUAL_STRING("(multiple)", row.basecaller_model);
    TEST_ASSERT_EQUAL_STRING("(multiple)", row.flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("(not meaningful)", row.run_date);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, -1.0, row.run_duration_s);
    TEST_ASSERT_EQUAL_UINT64(2, row.read_count);
    TEST_ASSERT_EQUAL_UINT64(300, row.base_count);

    fc_run_summary_destroy(model);
}

static void test_run_summary_uses_fallback_sample_and_earliest_time(void)
{
    fc_run_summary *model = fc_run_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_run_summary_read(
        model, 0, 100,
        "",
        "run-a",
        "model-a",
        "flow-a",
        "alias-a",
        "2026-01-01T02:00:00Z");
    add_run_summary_read(
        model, 0, 200,
        "",
        "run-a",
        "model-a",
        "flow-a",
        "alias-a",
        "2026-01-01T00:30:00Z");

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, NULL, "fallback-sample", "unclassified", &row));

    TEST_ASSERT_EQUAL_STRING("fallback-sample", row.sample_id);
    TEST_ASSERT_EQUAL_STRING("unclassified", row.barcode_id);
    TEST_ASSERT_EQUAL_STRING("alias-a", row.barcode_alias);
    TEST_ASSERT_EQUAL_STRING("run-a", row.run_id);
    TEST_ASSERT_EQUAL_STRING("model-a", row.basecaller_model);
    TEST_ASSERT_EQUAL_STRING("flow-a", row.flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("2026-01-01 00:30:00", row.run_date);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 5400.0, row.run_duration_s);
    TEST_ASSERT_EQUAL_UINT64(2, row.read_count);
    TEST_ASSERT_EQUAL_UINT64(300, row.base_count);

    fc_run_summary_destroy(model);
}

void run_summary_suite(void)
{
    RUN_TEST(test_run_summary_collapses_multiple_values);
    RUN_TEST(test_run_summary_uses_fallback_sample_and_earliest_time);
}
