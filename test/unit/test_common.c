#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "unity.h"
#include "main.h"
#include "core/fc_common.h"

FC_TEST(common_suite, test_min_max)
{
    int a = 5, b = 3;
    TEST_ASSERT_EQUAL_INT(3, min(a, b));
    TEST_ASSERT_EQUAL_INT(5, max(a, b));
    TEST_ASSERT_EQUAL_INT(10, max(min(10, 20), 5));
}

FC_TEST(common_suite, test_replace_char)
{
    char buf[32] = "A-C-G-T";
    int n = replace_char(buf, '-', '+');
    TEST_ASSERT_EQUAL_INT(3, n);
    TEST_ASSERT_EQUAL_STRING("A+C+G+T", buf);
}

FC_TEST(common_suite, test_create_rg_info)
{
    char rg[] = "12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode03";
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0", info->basecaller);
    TEST_ASSERT_NULL(info->modcaller);
    TEST_ASSERT_EQUAL_STRING("barcode03", info->barcode);
    destroy_rg_info(info);
}

FC_TEST(common_suite, test_create_rg_info_extra_at_without_modcaller_separator)
{
    char rg[] = "12345678-1234-1234-1234-123456789012_sup@v4.0.0@modcaller";
    readgroup* info = create_rg_info(rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0@modcaller", info->basecaller);
    TEST_ASSERT_NULL(info->modcaller);
    TEST_ASSERT_NULL(info->barcode);
    destroy_rg_info(info);
}

FC_TEST(common_suite, test_kahan_sum)
{
    double sum = 0.0, c = 0.0;
    for (int i = 0; i < 100000; ++i) {
        kahan_sum(&sum, 0.1, &c);
    }
    TEST_ASSERT_FLOAT_WITHIN(0.05, sum, 10000.0);
}

FC_TEST(common_suite, test_kahan_sum_cancellation)
{
    double sum = 0.0, c = 0.0;
    kahan_sum(&sum, 1e12, &c);
    for (int i = 0; i < 10; ++i) {
        kahan_sum(&sum, 0.1, &c);
    }
    kahan_sum(&sum, -1e12, &c);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, sum);
}

FC_TEST(common_suite, test_parse_int_strict)
{
    int v = 0;
    TEST_ASSERT_TRUE(fc_parse_int_strict("42", &v));
    TEST_ASSERT_EQUAL_INT(42, v);
    TEST_ASSERT_FALSE(fc_parse_int_strict("42x", &v));
    TEST_ASSERT_FALSE(fc_parse_int_strict("", &v));
}

FC_TEST(common_suite, test_parse_size_strict)
{
    size_t v = 0;
    TEST_ASSERT_TRUE(fc_parse_size_strict("123", &v));
    TEST_ASSERT_EQUAL_UINT64(123, (uint64_t)v);
    TEST_ASSERT_FALSE(fc_parse_size_strict("-1", &v));
    TEST_ASSERT_FALSE(fc_parse_size_strict("1.5", &v));
}

FC_TEST(common_suite, test_parse_u32_strict_overflow)
{
    uint32_t v = 0;
    TEST_ASSERT_TRUE(fc_parse_u32_strict("4294967295", &v));
    TEST_ASSERT_EQUAL_UINT32(4294967295u, v);
    TEST_ASSERT_FALSE(fc_parse_u32_strict("4294967296", &v));
}

FC_TEST(common_suite, test_parse_float_strict)
{
    float v = 0.0f;
    TEST_ASSERT_TRUE(fc_parse_float_strict("3.5", &v));
    TEST_ASSERT_FLOAT_WITHIN(0.0001f, 3.5f, v);
    TEST_ASSERT_FALSE(fc_parse_float_strict("nan", &v));
    TEST_ASSERT_FALSE(fc_parse_float_strict("inf", &v));
    TEST_ASSERT_FALSE(fc_parse_float_strict("3.5x", &v));
}
