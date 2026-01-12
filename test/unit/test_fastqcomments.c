#include <stdio.h>
#include <string.h>

#include "unity.h"
#include "fastqcomments.h"
#include "htslib/kstring.h"

static void test_parse_simple_comment(void) {
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "runid=abc123 barcode=barcode12 start_time=2023-01-01T00:00:00Z");

    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_EQUAL_STRING("abc123", meta->runid);
    TEST_ASSERT_EQUAL_STRING("barcode12", meta->barcode);
    TEST_ASSERT_EQUAL_INT(12, meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", meta->start_time);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RD:Z:abc123"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "BC:Z:barcode12"));

    destroy_read_meta(meta);
    free(ks.s);
}

static void test_original_comment_preserved(void) {
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "raw_comment\twith\tspecials");
    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_NOT_NULL(meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw_comment\twith\tspecials", meta->orig_comment);
    destroy_read_meta(meta);
    free(ks.s);
}

static void test_unknown_key_goes_to_co(void) {
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "runid=abc123 unknown=foo");

    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RD:Z:abc123"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "CO:Z:unknown=foo"));
    destroy_read_meta(meta);
    free(ks.s);
}

static void test_sam_style_tabs(void) {
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "RG:Z:rg1\tBC:Z:barcode01\tZZ:Z:extra");
    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_EQUAL_STRING("rg1", meta->rg);
    TEST_ASSERT_EQUAL_STRING("barcode01", meta->barcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:rg1"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "BC:Z:barcode01"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "ZZ:Z:extra"));
    destroy_read_meta(meta);
    free(ks.s);
}

static void test_barcode_unclassified(void) {
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "barcode=unclassified");

    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_EQUAL_STRING("unclassified", meta->barcode);
    TEST_ASSERT_EQUAL_INT(0, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "BC:Z:unclassified"));

    destroy_read_meta(meta);
    free(ks.s);
}

static void test_rg_populates_runid_and_basecaller(void) {
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v1";
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "RG:Z:%s", rg);

    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", meta->basecaller);
    TEST_ASSERT_NULL(meta->rg_info->barcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RD:Z:12345678-1234-1234-1234-123456789012"));

    destroy_read_meta(meta);
    free(ks.s);
}

static void test_rg_fallback(void) {
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode07";
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "RG:Z:%s\tBC:Z:barcode07", rg);

    read_meta meta = parse_read_meta(ks);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0", meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RD:Z:12345678-1234-1234-1234-123456789012"));

    destroy_read_meta(meta);
    free(ks.s);
}

void fastqcomments_suite(void) {
    RUN_TEST(test_parse_simple_comment);
    RUN_TEST(test_original_comment_preserved);
    RUN_TEST(test_unknown_key_goes_to_co);
    RUN_TEST(test_sam_style_tabs);
    RUN_TEST(test_barcode_unclassified);
    RUN_TEST(test_rg_populates_runid_and_basecaller);
    RUN_TEST(test_rg_fallback);
}
