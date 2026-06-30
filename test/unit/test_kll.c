#include <math.h>
#include <stdio.h>

#include "unity.h"
#include "core/fc_kll.h"

static void test_empty_sketch(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    TEST_ASSERT_EQUAL_size_t(0, fc_kll_count(sketch));

    double q50 = fc_kll_compute_quantile(sketch, 0.5);
    TEST_ASSERT_TRUE(isnan(q50));

    fc_kll_destroy(sketch);
}

static void test_single_value(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    fc_kll_add(sketch, 42.0);

    TEST_ASSERT_EQUAL_size_t(1, fc_kll_count(sketch));

    double q50 = fc_kll_compute_quantile(sketch, 0.5);
    TEST_ASSERT_FLOAT_WITHIN(0.1, 42.0, q50);

    fc_kll_destroy(sketch);
}

static void test_bad_quantile_input(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    // add some data
    for (int i = 0; i < 10; i++) {
        fc_kll_add(sketch, (double)i);
    }

    // test out-of-range quantile values
    double bad_low = fc_kll_compute_quantile(sketch, -0.1);
    TEST_ASSERT_TRUE(isnan(bad_low));

    double bad_high = fc_kll_compute_quantile(sketch, 1.5);
    TEST_ASSERT_TRUE(isnan(bad_high));

    // valid edge cases should work
    double q0 = fc_kll_compute_quantile(sketch, 0.0);
    double q1 = fc_kll_compute_quantile(sketch, 1.0);
    TEST_ASSERT_FALSE(isnan(q0));
    TEST_ASSERT_FALSE(isnan(q1));

    fc_kll_destroy(sketch);
}

static void test_nan_input(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    // NaN should be ignored
    fc_kll_add(sketch, NAN);
    TEST_ASSERT_EQUAL_size_t(0, fc_kll_count(sketch));

    // add valid value, then NaN
    fc_kll_add(sketch, 42.0);
    fc_kll_add(sketch, NAN);
    TEST_ASSERT_EQUAL_size_t(1, fc_kll_count(sketch));

    double q50 = fc_kll_compute_quantile(sketch, 0.5);
    TEST_ASSERT_FLOAT_WITHIN(0.1, 42.0, q50);

    fc_kll_destroy(sketch);
}

static void test_infinity_input(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    // infinity should be ignored (non-finite)
    fc_kll_add(sketch, INFINITY);
    fc_kll_add(sketch, -INFINITY);
    TEST_ASSERT_EQUAL_size_t(0, fc_kll_count(sketch));

    // mixed with valid values
    fc_kll_add(sketch, 10.0);
    fc_kll_add(sketch, INFINITY);
    fc_kll_add(sketch, 20.0);
    TEST_ASSERT_EQUAL_size_t(2, fc_kll_count(sketch));

    double q50 = fc_kll_compute_quantile(sketch, 0.5);
    TEST_ASSERT_FLOAT_WITHIN(5.0, 15.0, q50);

    fc_kll_destroy(sketch);
}

static void test_quantile_triple(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    // add 100 sequential values
    for (int i = 0; i < 100; i++) {
        fc_kll_add(sketch, (double)i);
    }

    fc_kll_quantile_triple triple = fc_kll_compute_quantile_triple(sketch);

    // verify ordering
    TEST_ASSERT_TRUE(triple.q25 <= triple.q50);
    TEST_ASSERT_TRUE(triple.q50 <= triple.q75);

    // verify rough correctness (exact mode with 100 items)
    TEST_ASSERT_FLOAT_WITHIN(2.0, 25.0, triple.q25);
    TEST_ASSERT_FLOAT_WITHIN(2.0, 50.0, triple.q50);
    TEST_ASSERT_FLOAT_WITHIN(2.0, 75.0, triple.q75);

    fc_kll_destroy(sketch);
}

static void test_quantile_triple_empty(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    fc_kll_quantile_triple triple = fc_kll_compute_quantile_triple(sketch);

    // empty sketch should return NaN
    TEST_ASSERT_TRUE(isnan(triple.q25));
    TEST_ASSERT_TRUE(isnan(triple.q50));
    TEST_ASSERT_TRUE(isnan(triple.q75));

    fc_kll_destroy(sketch);
}

void kll_suite(void)
{
    RUN_TEST(test_empty_sketch);
    RUN_TEST(test_single_value);
    RUN_TEST(test_bad_quantile_input);
    RUN_TEST(test_nan_input);
    RUN_TEST(test_infinity_input);
    RUN_TEST(test_quantile_triple);
    RUN_TEST(test_quantile_triple_empty);
}
