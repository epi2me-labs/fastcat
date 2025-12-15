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

static void test_substring(void) {
    char src[] = "abcdefgh";
    char* sub = substring(src, 2, 3);
    TEST_ASSERT_EQUAL_STRING("cde", sub);
    free(sub);
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

void common_suite(void) {
    RUN_TEST(test_min_max);
    RUN_TEST(test_replace_char);
    RUN_TEST(test_substring);
    RUN_TEST(test_create_rg_info);
    RUN_TEST(test_kahan_sum);
}
