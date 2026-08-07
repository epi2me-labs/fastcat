#include <string.h>

#include "unity.h"
#include "main.h"

#include "htslib/sam.h"

#include "core/fc_bam_rg_cache.h"

typedef struct {
    const char *id;
    const char *ds;
    const char *runid;
    const char *basecall_model;
    const char *modbase_models;
    const char *dt;
    const char *lb;
    const char *pu;
    const char *pm;
    const char *sm;
    const char *al;
} rg_expectation;

typedef struct {
    const char *header;
    const rg_expectation *entries;
    size_t n_entries;
    const char *missing_id;
} rg_cache_fixture;

static void assert_optional_string(const char *expected, const char *actual)
{
    if (NULL == expected) {
        TEST_ASSERT_NULL(actual);
        return;
    }
    TEST_ASSERT_EQUAL_STRING(expected, actual);
}

static void assert_rg_entry(const fc_bam_rg_entry *entry, const rg_expectation *expected)
{
    TEST_ASSERT_NOT_NULL(entry);
    TEST_ASSERT_EQUAL_STRING(expected->id, entry->id);
    assert_optional_string(expected->ds, entry->ds);
    assert_optional_string(expected->runid, entry->runid);
    assert_optional_string(expected->basecall_model, entry->basecall_model);
    assert_optional_string(expected->modbase_models, entry->modbase_models);
    assert_optional_string(expected->dt, entry->dt);
    assert_optional_string(expected->lb, entry->lb);
    assert_optional_string(expected->pu, entry->pu);
    assert_optional_string(expected->pm, entry->pm);
    assert_optional_string(expected->sm, entry->sm);
    assert_optional_string(expected->al, entry->al);
}

static void run_rg_cache_fixture(const rg_cache_fixture *fixture)
{
    bam_hdr_t *hdr = sam_hdr_parse((int)strlen(fixture->header), fixture->header);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_rg_cache *cache = fc_bam_rg_cache_create(hdr);
    TEST_ASSERT_NOT_NULL(cache);

    for (size_t i = 0; i < fixture->n_entries; ++i) {
        const rg_expectation *expected = &fixture->entries[i];
        assert_rg_entry(fc_bam_rg_cache_find(cache, expected->id), expected);
    }

    if (NULL != fixture->missing_id) {
        TEST_ASSERT_NULL(fc_bam_rg_cache_find(cache, fixture->missing_id));
    }

    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_rg_cache_suite, test_bam_rg_cache_parses_basecall_model_only)
{
    static const rg_expectation entries[] = {
        {
            .id = "rg1",
            .ds = "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model=rna004_130bps_fast@v5.1.0",
            .runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5",
            .basecall_model = "rna004_130bps_fast@v5.1.0",
        },
    };
    static const rg_cache_fixture fixture = {
        .header =
            "@HD\tVN:1.6\n"
            "@RG\tID:rg1\tDS:runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model=rna004_130bps_fast@v5.1.0\n",
        .entries = entries,
        .n_entries = sizeof(entries) / sizeof(entries[0]),
    };

    run_rg_cache_fixture(&fixture);
}

FC_TEST(bam_rg_cache_suite, test_bam_rg_cache_parses_modbase_models_variants)
{
    static const rg_expectation entries[] = {
        {
            .id = "rg1",
            .ds = "basecall_model=rna004_130bps_sup@2025-01-22_n1mpU modbase_models=modbase_model_version_id",
            .basecall_model = "rna004_130bps_sup@2025-01-22_n1mpU",
            .modbase_models = "modbase_model_version_id",
        },
        {
            .id = "rg2",
            .ds = "basecall_model=rna004_130bps_fast@v5.1.0 modbase_models=rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1",
            .basecall_model = "rna004_130bps_fast@v5.1.0",
            .modbase_models = "rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1",
        },
    };
    static const rg_cache_fixture fixture = {
        .header =
            "@HD\tVN:1.6\n"
            "@RG\tID:rg1\tDS:basecall_model=rna004_130bps_sup@2025-01-22_n1mpU modbase_models=modbase_model_version_id\n"
            "@RG\tID:rg2\tDS:basecall_model=rna004_130bps_fast@v5.1.0 modbase_models=rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1\n",
        .entries = entries,
        .n_entries = sizeof(entries) / sizeof(entries[0]),
    };

    run_rg_cache_fixture(&fixture);
}

FC_TEST(bam_rg_cache_suite, test_bam_rg_cache_lookups_and_rg_fields)
{
    static const rg_expectation entries[] = {
        {
            .id = "rg1",
            .ds = "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model=rna004_130bps_sup@2024-06-05_canonical modbase_models=modbase_model_version_id",
            .runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5",
            .basecall_model = "rna004_130bps_sup@2024-06-05_canonical",
            .modbase_models = "modbase_model_version_id",
            .dt = "2025-01-13T10:45:28.681306+00:00",
            .lb = "sampleA",
            .pu = "FLO-MIN114",
            .pm = "MIN-001",
            .sm = "barcode07",
            .al = "alias07",
        },
        {
            .id = "rg2",
            .ds = "basecall_model=rna004_130bps_sup@2025-06-12_n1mpU-integrity",
            .basecall_model = "rna004_130bps_sup@2025-06-12_n1mpU-integrity",
            .dt = "2016-01-19T15:21:32.59+02:00",
            .lb = "sampleB",
            .pu = "FBH03097",
            .pm = "PRO-002",
            .sm = "unclassified",
            .al = "unknown",
        },
    };
    static const rg_cache_fixture fixture = {
        .header =
            "@HD\tVN:1.6\n"
            "@RG\tID:rg1\tDT:2025-01-13T10:45:28.681306+00:00\tLB:sampleA\tPU:FLO-MIN114\tPM:MIN-001\tSM:barcode07\tal:alias07\tDS:runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model=rna004_130bps_sup@2024-06-05_canonical modbase_models=modbase_model_version_id\n"
            "@RG\tID:rg2\tDT:2016-01-19T15:21:32.59+02:00\tLB:sampleB\tPU:FBH03097\tPM:PRO-002\tSM:unclassified\tal:unknown\tDS:basecall_model=rna004_130bps_sup@2025-06-12_n1mpU-integrity\n",
        .entries = entries,
        .n_entries = sizeof(entries) / sizeof(entries[0]),
        .missing_id = "missing",
    };

    run_rg_cache_fixture(&fixture);
}

FC_TEST(bam_rg_cache_suite, test_bam_rg_cache_handles_partial_or_missing_ds)
{
    static const rg_expectation entries[] = {
        {
            .id = "rg1",
            .ds = "runid=1234",
            .runid = "1234",
        },
        {
            .id = "rg2",
            .ds = "modbase_models=foo@v1,bar@v2",
            .modbase_models = "foo@v1,bar@v2",
        },
        {
            .id = "rg3",
            .lb = "sampleC",
        },
    };
    static const rg_cache_fixture fixture = {
        .header =
            "@HD\tVN:1.6\n"
            "@RG\tID:rg1\tDS:runid=1234\n"
            "@RG\tID:rg2\tDS:modbase_models=foo@v1,bar@v2\n"
            "@RG\tID:rg3\tLB:sampleC\n",
        .entries = entries,
        .n_entries = sizeof(entries) / sizeof(entries[0]),
    };

    run_rg_cache_fixture(&fixture);
}

FC_TEST(bam_rg_cache_suite, test_bam_rg_cache_reuses_legacy_rg_parsing)
{
    const char *header =
        "@HD\tVN:1.6\n"
        "@RG\tID:e4994c62-93f9-439a-bc8f-d20c95a137a5_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n";
    bam_hdr_t *hdr = sam_hdr_parse((int)strlen(header), header);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_rg_cache *cache = fc_bam_rg_cache_create(hdr);
    TEST_ASSERT_NOT_NULL(cache);

    const char *header_rg = "e4994c62-93f9-439a-bc8f-d20c95a137a5_dna_r10.4.1_e8.2_400bps_sup@v5.0.0";
    const fc_bam_rg_entry *first = fc_bam_rg_cache_resolve(cache, header_rg);
    const fc_bam_rg_entry *second = fc_bam_rg_cache_resolve(cache, header_rg);
    TEST_ASSERT_NOT_NULL(first);
    TEST_ASSERT_EQUAL_PTR(first, second);
    TEST_ASSERT_NOT_NULL(first->parsed_id);
    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5", first->parsed_id->runid);

    const char *unknown_rg = "d7dbdc6e-f26b-4d7f-995f-b7034bcf7290_dna_r10.4.1_e8.2_400bps_hac@v5.0.0";
    const fc_bam_rg_entry *unknown_first = fc_bam_rg_cache_resolve(cache, unknown_rg);
    const fc_bam_rg_entry *unknown_second = fc_bam_rg_cache_resolve(cache, unknown_rg);
    TEST_ASSERT_NOT_NULL(unknown_first);
    TEST_ASSERT_EQUAL_PTR(unknown_first, unknown_second);
    TEST_ASSERT_NOT_NULL(unknown_first->parsed_id);

    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}
