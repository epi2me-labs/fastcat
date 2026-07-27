#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "unity.h"
#include "main.h"

#include "core/fc_run_summary.h"
#include "core/fc_stats_engine.h"

static void add_run_summary_read(
    fc_run_summary *model,
    size_t barcode,
    size_t length,
    const char *sample_id,
    const char *run_id,
    const char *basecaller_model,
    const char *modbase_models,
    const char *flow_cell_id,
    const char *barcode_alias,
    const char *start_time)
{
    bam1_t bam = { 0 };
    const char *meta_sample_id = "";
    const char *meta_run_id = "";
    const char *meta_basecaller_model = "";
    const char *meta_modbase_models = "";
    const char *meta_flow_cell_id = "";
    const char *meta_barcode_alias = "";
    const char *meta_start_time = "";

    if (NULL != sample_id) {
        meta_sample_id = sample_id;
    }
    if (NULL != run_id) {
        meta_run_id = run_id;
    }
    if (NULL != basecaller_model) {
        meta_basecaller_model = basecaller_model;
    }
    if (NULL != modbase_models) {
        meta_modbase_models = modbase_models;
    }
    if (NULL != flow_cell_id) {
        meta_flow_cell_id = flow_cell_id;
    }
    if (NULL != barcode_alias) {
        meta_barcode_alias = barcode_alias;
    }
    if (NULL != start_time) {
        meta_start_time = start_time;
    }

    read_meta meta = {
        .sample_id = (char *)meta_sample_id,
        .runid = (char *)meta_run_id,
        .basecaller = (char *)meta_basecaller_model,
        .modcallers = (char *)meta_modbase_models,
        .flow_cell_id = (char *)meta_flow_cell_id,
        .barcode_alias = (char *)meta_barcode_alias,
        .start_time = (char *)meta_start_time
    };
    fc_read read = {
        .bam = &bam,
        .meta = &meta
    };

    bam.core.l_qseq = (int32_t)length;
    fc_run_summary_add_read(model, barcode, &read);
}

static void assert_joined_pair(const char *actual, const char *a, const char *b)
{
    // we join per-read values in arbitrary order, so either order is valid.
    // They are joined with a semi-colon because the output spec already joins
    // multiple values per-read with a comma. Our output is the unique set of
    // strings seen, not the unique set of models (possibly different order).
    char expected_ab[1024];
    char expected_ba[1024];
    snprintf(expected_ab, sizeof(expected_ab), "%s;%s", a, b);
    snprintf(expected_ba, sizeof(expected_ba), "%s;%s", b, a);
    TEST_ASSERT_TRUE(0 == strcmp(actual, expected_ab) || 0 == strcmp(actual, expected_ba));
}

#define MODBASE_MODELS_1 "rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1"
#define MODBASE_MODELS_2 "rna004_130bps_hac@v5.3.0_m5C@v1,rna004_130bps_hac@v5.4.0_pseU@v2"

FC_TEST(run_summary_suite, test_run_summary_collapses_multiple_values)
{
    fc_run_summary *model = fc_run_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_run_summary_read(
        model, 0, 100,
        "sample-a",
        "run-a",
        "model-a",
        MODBASE_MODELS_1,
        "flow-a",
        "alias-a",
        "2026-01-01T00:00:00Z");
    add_run_summary_read(
        model, 0, 200,
        "sample-b",
        "run-b",
        "model-b",
        MODBASE_MODELS_2,
        "flow-b",
        "alias-b",
        "2026-01-01T01:00:00Z");

    TEST_ASSERT_TRUE(fc_run_summary_seen(model, 0));

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, NULL, false, NULL, "fallback-sample", "barcode0001", &row));

    assert_joined_pair(row.sample_id, "sample-a", "sample-b");
    TEST_ASSERT_EQUAL_STRING("fastq", row.input_type);
    TEST_ASSERT_EQUAL_STRING("barcode0001", row.barcode_id);
    assert_joined_pair(row.barcode_alias, "alias-a", "alias-b");
    assert_joined_pair(row.run_id, "run-a", "run-b");
    assert_joined_pair(row.basecaller_model, "model-a", "model-b");
    assert_joined_pair(row.modbase_models, MODBASE_MODELS_1, MODBASE_MODELS_2);
    assert_joined_pair(row.flow_cell_id, "flow-a", "flow-b");
    TEST_ASSERT_EQUAL_STRING("(not meaningful)", row.run_date);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, -1.0, row.run_duration_s);
    TEST_ASSERT_EQUAL_UINT64(2, row.read_count);
    TEST_ASSERT_EQUAL_UINT64(300, row.base_count);

    fc_run_summary_row_destroy(&row);
    fc_run_summary_destroy(model);
}

FC_TEST(run_summary_suite, test_run_summary_uses_fallback_sample_and_earliest_time)
{
    fc_run_summary *model = fc_run_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_run_summary_read(
        model, 0, 100,
        "",
        "run-a",
        "model-a",
        "mod-a",
        "flow-a",
        "alias-a",
        "2026-01-01T02:00:00Z");
    add_run_summary_read(
        model, 0, 200,
        "",
        "run-a",
        "model-a",
        "mod-a",
        "flow-a",
        "alias-a",
        "2026-01-01T00:30:00Z");

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, NULL, false, NULL, "fallback-sample", "unclassified", &row));

    TEST_ASSERT_EQUAL_STRING("fallback-sample", row.sample_id);
    TEST_ASSERT_EQUAL_STRING("fastq", row.input_type);
    TEST_ASSERT_EQUAL_STRING("unclassified", row.barcode_id);
    TEST_ASSERT_EQUAL_STRING("alias-a", row.barcode_alias);
    TEST_ASSERT_EQUAL_STRING("run-a", row.run_id);
    TEST_ASSERT_EQUAL_STRING("model-a", row.basecaller_model);
    TEST_ASSERT_EQUAL_STRING("mod-a", row.modbase_models);
    TEST_ASSERT_EQUAL_STRING("flow-a", row.flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("2026-01-01 00:30:00", row.run_date);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 5400.0, row.run_duration_s);
    TEST_ASSERT_EQUAL_UINT64(2, row.read_count);
    TEST_ASSERT_EQUAL_UINT64(300, row.base_count);

    fc_run_summary_row_destroy(&row);
    fc_run_summary_destroy(model);
}

FC_TEST(run_summary_suite, test_run_summary_classifies_bam_input_type)
{
    fc_run_summary *model = fc_run_summary_create(1);
    fc_stats_engine *engine = fc_stats_engine_create(false, false, 0.0f, 0.0f, false);
    TEST_ASSERT_NOT_NULL(model);
    TEST_ASSERT_NOT_NULL(engine);

    add_run_summary_read(
        model, 0, 100,
        "sample-a",
        "run-a",
        "model-a",
        "",
        "flow-a",
        "alias-a",
        "2026-01-01T00:00:00Z");

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, engine, true, NULL, "fallback-sample", "unclassified", &row));
    TEST_ASSERT_EQUAL_STRING("unaligned bam", row.input_type);
    fc_run_summary_row_destroy(&row);

    add_length_count(engine->length_stats, 100);
    TEST_ASSERT_EQUAL_INT(
        0,
        fc_run_summary_materialize(model, 0, engine, true, NULL, "fallback-sample", "unclassified", &row));
    TEST_ASSERT_EQUAL_STRING("aligned bam", row.input_type);

    fc_run_summary_row_destroy(&row);
    fc_stats_engine_destroy(engine);
    fc_run_summary_destroy(model);
}

FC_TEST(run_summary_suite, test_run_summary_materializes_filter_counts_per_barcode)
{
    fc_run_summary *model = fc_run_summary_create(2);
    TEST_ASSERT_NOT_NULL(model);

    fc_run_summary_add_filter_event(model, 1, R_RECORD_SEEN);
    fc_run_summary_add_filter_event(model, 1, R_RECORD_SEEN);
    fc_run_summary_add_filter_event(model, 1, R_RECORD_SEEN);
    fc_run_summary_add_filter_event(model, 1, R_RECORD_OK);
    fc_run_summary_add_filter_event(model, 1, R_TOO_SHORT);
    fc_run_summary_add_filter_event(model, 1, R_LOW_QUALITY);
    fc_run_summary_add_filter_event(model, 0, R_RECORD_SEEN);

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(0,
        fc_run_summary_materialize(model, 1, NULL, false, NULL, "sample", "barcode0001", &row));
    TEST_ASSERT_EQUAL_UINT64(3, row.records_seen);
    TEST_ASSERT_EQUAL_UINT64(1, row.records_passed_filters);
    TEST_ASSERT_EQUAL_UINT64(0, row.records_too_long);
    TEST_ASSERT_EQUAL_UINT64(1, row.records_too_short);
    TEST_ASSERT_EQUAL_UINT64(1, row.records_low_quality);
    TEST_ASSERT_EQUAL_UINT64(0, row.records_dust_masked);

    fc_run_summary_row_destroy(&row);
    fc_run_summary_destroy(model);
}

FC_TEST(run_summary_suite, test_run_summary_nxx_lxx_metrics)
{
    fc_run_summary *model = fc_run_summary_create(1);
    read_stats *length_stats = create_length_stats();
    TEST_ASSERT_NOT_NULL(model);
    TEST_ASSERT_NOT_NULL(length_stats);

    add_run_summary_read(model, 0, 100, "sample", "run", "model", "", "flow", "alias", "");
    add_run_summary_read(model, 0, 200, "sample", "run", "model", "", "flow", "alias", "");
    add_run_summary_read(model, 0, 400, "sample", "run", "model", "", "flow", "alias", "");
    add_length_count(length_stats, 100);
    add_length_count(length_stats, 200);
    add_length_count(length_stats, 400);

    fc_run_summary_row row = { 0 };
    TEST_ASSERT_EQUAL_INT(0,
        fc_run_summary_materialize(model, 0, NULL, false, length_stats, "fallback", "all", &row));

    TEST_ASSERT_EQUAL_UINT64(400, row.n5_bp);
    TEST_ASSERT_EQUAL_UINT64(400, row.n10_bp);
    TEST_ASSERT_EQUAL_UINT64(400, row.n25_bp);
    TEST_ASSERT_EQUAL_UINT64(400, row.n50_bp);
    TEST_ASSERT_EQUAL_UINT64(200, row.n75_bp);
    TEST_ASSERT_EQUAL_UINT64(100, row.n90_bp);
    TEST_ASSERT_EQUAL_UINT64(100, row.n95_bp);
    TEST_ASSERT_EQUAL_UINT64(1, row.l5);
    TEST_ASSERT_EQUAL_UINT64(1, row.l10);
    TEST_ASSERT_EQUAL_UINT64(1, row.l25);
    TEST_ASSERT_EQUAL_UINT64(1, row.l50);
    TEST_ASSERT_EQUAL_UINT64(2, row.l75);
    TEST_ASSERT_EQUAL_UINT64(3, row.l90);
    TEST_ASSERT_EQUAL_UINT64(3, row.l95);

    fc_run_summary_row_destroy(&row);
    destroy_length_stats(length_stats);
    fc_run_summary_destroy(model);
}
