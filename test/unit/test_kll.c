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

static void test_quantile_linear_interpolation(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);
    fc_kll_add(sketch, 0.0);
    fc_kll_add(sketch, 10.0);

    double q50 = fc_kll_compute_quantile(sketch, 0.5);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 5.0, q50);

    double ps[2] = { 0.25, 0.75 };
    double out[2] = { NAN, NAN };
    int rc = fc_kll_compute_quantiles(sketch, ps, 2, out);

    TEST_ASSERT_EQUAL_INT(0, rc);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2.5, out[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 7.5, out[1]);

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

static void test_density_from_quantiles_basic(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);
    fc_kll_add(sketch, 0.0);
    fc_kll_add(sketch, 1.0);
    fc_kll_add(sketch, 2.0);

    fc_kll_density density = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 3, &density, false, 1));

    TEST_ASSERT_EQUAL_size_t(2, density.n);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 0.5f, (float)density.support[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 1.5f, (float)density.support[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 0.5f, (float)density.values[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 0.5f, (float)density.values[1]);

    // For p-grid [0, 0.5, 1.0], interval widths are both 1.0.
    double integral = density.values[0] + density.values[1];
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 1.0f, (float)integral);

    fc_kll_density_destroy(&density);
    fc_kll_destroy(sketch);
}

static void test_density_from_quantiles_smooth_changes_values(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);

    fc_kll_add(sketch, 0.0);
    for (int i = 0; i < 8; i++) {
        fc_kll_add(sketch, 1.0);
    }
    fc_kll_add(sketch, 10.0);

    fc_kll_density raw = { 0 };
    fc_kll_density smooth = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 33, &raw, false, 1));
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 33, &smooth, true, 3));

    TEST_ASSERT_EQUAL_size_t(raw.n, smooth.n);
    TEST_ASSERT_TRUE(raw.n > 0);

    int any_diff = 0;
    for (size_t i = 0; i < raw.n; i++) {
        TEST_ASSERT_FLOAT_WITHIN(1e-9f, (float)raw.support[i], (float)smooth.support[i]);
        TEST_ASSERT_TRUE(isfinite(smooth.values[i]));
        TEST_ASSERT_TRUE(smooth.values[i] >= 0.0);
        if (fabs(smooth.values[i] - raw.values[i]) > 1e-12) {
            any_diff = 1;
        }
    }
    TEST_ASSERT_TRUE(any_diff);

    fc_kll_density_destroy(&raw);
    fc_kll_density_destroy(&smooth);
    fc_kll_destroy(sketch);
}

static void test_density_from_quantiles_all_equal(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);
    for (int i = 0; i < 10; i++) {
        fc_kll_add(sketch, 7.0);
    }

    fc_kll_density density = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 17, &density, true, 3));

    TEST_ASSERT_EQUAL_size_t(1, density.n);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 7.0f, (float)density.support[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-9f, 1.0f, (float)density.values[0]);

    fc_kll_density_destroy(&density);
    fc_kll_destroy(sketch);
}

static void test_density_from_quantiles_smoothing_window_size(void)
{
    fc_kll_sketch *sketch = fc_kll_create(200);
    fc_kll_add(sketch, 0.0);
    for (int i = 0; i < 12; i++) {
        fc_kll_add(sketch, 1.0);
    }
    for (int i = 0; i < 6; i++) {
        fc_kll_add(sketch, 2.0);
    }
    fc_kll_add(sketch, 10.0);

    fc_kll_density smooth3 = { 0 };
    fc_kll_density smooth9 = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 65, &smooth3, true, 3));
    TEST_ASSERT_EQUAL_INT(0, fc_kll_compute_density_from_quantiles(sketch, 65, &smooth9, true, 9));

    TEST_ASSERT_EQUAL_size_t(smooth3.n, smooth9.n);
    TEST_ASSERT_TRUE(smooth3.n > 0);

    int any_diff = 0;
    for (size_t i = 0; i < smooth3.n; i++) {
        TEST_ASSERT_FLOAT_WITHIN(1e-9f, (float)smooth3.support[i], (float)smooth9.support[i]);
        if (fabs(smooth3.values[i] - smooth9.values[i]) > 1e-12) {
            any_diff = 1;
        }
    }
    TEST_ASSERT_TRUE(any_diff);

    // even smoothing window is invalid when smooth=true
    TEST_ASSERT_EQUAL_INT(-1, fc_kll_compute_density_from_quantiles(sketch, 65, &smooth3, true, 4));

    fc_kll_density_destroy(&smooth3);
    fc_kll_density_destroy(&smooth9);
    fc_kll_destroy(sketch);
}

void kll_suite(void)
{
    RUN_TEST(test_empty_sketch);
    RUN_TEST(test_single_value);
    RUN_TEST(test_bad_quantile_input);
    RUN_TEST(test_quantile_linear_interpolation);
    RUN_TEST(test_nan_input);
    RUN_TEST(test_infinity_input);
    RUN_TEST(test_quantile_triple);
    RUN_TEST(test_quantile_triple_empty);
    RUN_TEST(test_density_from_quantiles_basic);
    RUN_TEST(test_density_from_quantiles_smooth_changes_values);
    RUN_TEST(test_density_from_quantiles_all_equal);
    RUN_TEST(test_density_from_quantiles_smoothing_window_size);
}
