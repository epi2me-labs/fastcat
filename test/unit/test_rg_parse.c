#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "main.h"
#include "core/fc_common.h"

typedef struct {
    const char* runid;
    const char* basecall_model;
    const char* mod_model;
    const char* barcode;
    const char* suffix;
} rg_case;

static char* build_rg(const rg_case* c) {
    size_t len = 0;
    len += strlen(c->runid);
    if (c->basecall_model) len += 1 + strlen(c->basecall_model);
    if (c->mod_model) len += 1 + strlen(c->mod_model);
    if (c->barcode) len += 1 + strlen(c->barcode);
    if (c->suffix) len += strlen(c->suffix);

    char* rg = calloc(len + 1, 1);
    strcat(rg, c->runid);
    if (c->basecall_model) {
        strcat(rg, "_");
        strcat(rg, c->basecall_model);
    }
    if (c->mod_model) {
        strcat(rg, "_");
        strcat(rg, c->mod_model);
    }
    if (c->barcode) {
        strcat(rg, "_");
        strcat(rg, c->barcode);
    }
    if (c->suffix) {
        strcat(rg, c->suffix);
    }
    return rg;
}

FC_TEST(rg_parse_suite, test_rg_matrix) {
    const char *runid_acq = "ef1af1ab8967cb20ca30dbeca93fd66592bf4619";
    const char *runid_proto = "c886531d-28f5-41f6-b948-948e8cb78e5e";
    const char *basecall = "basecall_model_name@v1.2.3";
    const char *mod_model = "basecall_model_name@v1.2.3_5mCG_5hmCG@v1";
    const char *barcode = "barcode01";
    const char *suffix = "-1A2B3C4D";

    const rg_case cases[] = {
        {runid_acq, basecall, mod_model, barcode, suffix},
        {runid_acq, basecall, mod_model, barcode, NULL},
        {runid_acq, basecall, mod_model, NULL, suffix},
        {runid_acq, basecall, mod_model, NULL, NULL},
        {runid_acq, basecall, NULL, barcode, suffix},
        {runid_acq, basecall, NULL, barcode, NULL},
        {runid_acq, basecall, NULL, NULL, suffix},
        {runid_acq, basecall, NULL, NULL, NULL},
        {runid_proto, basecall, mod_model, barcode, suffix},
        {runid_proto, basecall, mod_model, barcode, NULL},
        {runid_proto, basecall, mod_model, NULL, suffix},
        {runid_proto, basecall, mod_model, NULL, NULL},
        {runid_proto, basecall, NULL, barcode, suffix},
        {runid_proto, basecall, NULL, barcode, NULL},
        {runid_proto, basecall, NULL, NULL, suffix},
        {runid_proto, basecall, NULL, NULL, NULL},
    };

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        char* rg = build_rg(&cases[i]);
        readgroup* info = create_rg_info(rg);

        TEST_ASSERT_EQUAL_STRING(cases[i].runid, info->runid);
        if (cases[i].basecall_model == NULL) {
            TEST_ASSERT_NULL(info->basecaller);
        } else {
            TEST_ASSERT_EQUAL_STRING(cases[i].basecall_model, info->basecaller);
        }
        if (cases[i].mod_model == NULL) {
            TEST_ASSERT_NULL(info->modcaller);
        } else {
            TEST_ASSERT_EQUAL_STRING(cases[i].mod_model, info->modcaller);
        }
        if (cases[i].barcode == NULL) {
            TEST_ASSERT_NULL(info->barcode);
        } else {
            TEST_ASSERT_EQUAL_STRING(cases[i].barcode, info->barcode);
        }

        destroy_rg_info(info);
        free(rg);
    }
}

FC_TEST(rg_parse_suite, test_hex_suffix_stripped) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_sup@v1_barcode01-ABCDEF");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", info->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode01", info->barcode);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_hex_suffix_no_barcode) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_sup@v1-ABCDEF");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", info->basecaller);
    TEST_ASSERT_NULL(info->barcode);
    TEST_ASSERT_NULL(info->modcaller);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_hex_suffix_limits) {
    char *eight_digits = strdup("12345678-1234-1234-1234-123456789012_sup@v1_barcode01-1A2B3C4D");
    readgroup *info = create_rg_info(eight_digits);
    TEST_ASSERT_EQUAL_STRING("barcode01", info->barcode);
    destroy_rg_info(info);
    free(eight_digits);

    char *nine_digits = strdup("12345678-1234-1234-1234-123456789012_sup@v1_barcode01-1A2B3C4D5");
    info = create_rg_info(nine_digits);
    TEST_ASSERT_EQUAL_STRING("barcode01-1A2B3C4D5", info->barcode);
    destroy_rg_info(info);
    free(nine_digits);

    char *no_digits = strdup("12345678-1234-1234-1234-123456789012_sup@v1_barcode01-");
    info = create_rg_info(no_digits);
    TEST_ASSERT_EQUAL_STRING("barcode01", info->barcode);
    destroy_rg_info(info);
    free(no_digits);
}

FC_TEST(rg_parse_suite, test_no_at_in_basecaller) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_supv1_barcode01");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_NULL(info->basecaller);
    TEST_ASSERT_NULL(info->modcaller);
    TEST_ASSERT_NULL(info->barcode);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_invalid_runid_length) {
    char* rg = strdup("short_sup@v1_barcode01");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_NULL(info->runid);
    TEST_ASSERT_NULL(info->basecaller);
    TEST_ASSERT_NULL(info->barcode);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_modcaller_extraction) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_sup@v1_mod@v2_barcode01");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("sup@v1", info->basecaller);
    TEST_ASSERT_EQUAL_STRING("mod@v2", info->modcaller);
    TEST_ASSERT_EQUAL_STRING("barcode01", info->barcode);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_modcaller_no_barcode) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_sup@v1_mod@v2");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("sup@v1", info->basecaller);
    TEST_ASSERT_EQUAL_STRING("mod@v2", info->modcaller);
    TEST_ASSERT_NULL(info->barcode);
    destroy_rg_info(info);
    free(rg);
}

FC_TEST(rg_parse_suite, test_no_barcode) {
    char* rg = strdup("12345678-1234-1234-1234-123456789012_sup@v1");
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("sup@v1", info->basecaller);
    TEST_ASSERT_NULL(info->barcode);
    destroy_rg_info(info);
    free(rg);
}
