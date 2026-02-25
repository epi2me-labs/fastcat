#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "unity.h"
#include "common.h"

static void test_min_max(void) {
    int a = 5, b = 3;
    TEST_ASSERT_EQUAL_INT(3, min(a, b));
    TEST_ASSERT_EQUAL_INT(5, max(a, b));
    TEST_ASSERT_EQUAL_INT(10, max(min(10, 20), 5));
}

static void test_replace_char(void) {
    char buf[32] = "A-C-G-T";
    int n = replace_char(buf, '-', '+');
    TEST_ASSERT_EQUAL_INT(3, n);
    TEST_ASSERT_EQUAL_STRING("A+C+G+T", buf);
}

static void test_create_rg_info(void) {
    char rg[] = "12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode03";
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0", info->basecaller);
    TEST_ASSERT_NULL(info->modcaller);
    TEST_ASSERT_EQUAL_STRING("barcode03", info->barcode);
    destroy_rg_info(info);
}

static void test_kahan_sum(void) {
    double sum = 0.0, c = 0.0;
    for (int i = 0; i < 100000; ++i) {
        kahan_sum(&sum, 0.1, &c);
    }
    TEST_ASSERT_FLOAT_WITHIN(0.05, sum, 10000.0);
}

static void test_parse_int_strict(void) {
    int v = 0;
    TEST_ASSERT_TRUE(fc_parse_int_strict("42", &v));
    TEST_ASSERT_EQUAL_INT(42, v);
    TEST_ASSERT_FALSE(fc_parse_int_strict("42x", &v));
    TEST_ASSERT_FALSE(fc_parse_int_strict("", &v));
}

static void test_parse_size_strict(void) {
    size_t v = 0;
    TEST_ASSERT_TRUE(fc_parse_size_strict("123", &v));
    TEST_ASSERT_EQUAL_UINT64(123, (uint64_t)v);
    TEST_ASSERT_FALSE(fc_parse_size_strict("-1", &v));
    TEST_ASSERT_FALSE(fc_parse_size_strict("1.5", &v));
}

static void test_parse_u32_strict_overflow(void) {
    uint32_t v = 0;
    TEST_ASSERT_TRUE(fc_parse_u32_strict("4294967295", &v));
    TEST_ASSERT_EQUAL_UINT32(4294967295u, v);
    TEST_ASSERT_FALSE(fc_parse_u32_strict("4294967296", &v));
}

static void test_parse_float_strict(void) {
    float v = 0.0f;
    TEST_ASSERT_TRUE(fc_parse_float_strict("3.5", &v));
    TEST_ASSERT_FLOAT_WITHIN(0.0001f, 3.5f, v);
    TEST_ASSERT_FALSE(fc_parse_float_strict("nan", &v));
    TEST_ASSERT_FALSE(fc_parse_float_strict("inf", &v));
    TEST_ASSERT_FALSE(fc_parse_float_strict("3.5x", &v));
}

void common_suite(void) {
    RUN_TEST(test_min_max);
    RUN_TEST(test_replace_char);
    RUN_TEST(test_create_rg_info);
    RUN_TEST(test_kahan_sum);
    RUN_TEST(test_parse_int_strict);
    RUN_TEST(test_parse_size_strict);
    RUN_TEST(test_parse_u32_strict_overflow);
    RUN_TEST(test_parse_float_strict);
}
