#include <math.h>
#include <stddef.h>

#include "unity.h"
#include "core/fc_length_qscore_2d.h"

static double matrix_sum(const fc_length_qscore_2d_matrix *matrix)
{
    double sum = 0.0;
    for (size_t i = 0; i < matrix->length_bins * matrix->qscore_bins; i++) {
        sum += matrix->counts[i];
    }
    return sum;
}

static double column_sum(const fc_length_qscore_2d_matrix *matrix, size_t qbin)
{
    double sum = 0.0;
    const double *col = matrix->counts + qbin * matrix->length_bins;
    for (size_t i = 0; i < matrix->length_bins; i++) {
        sum += col[i];
    }
    return sum;
}

static void assert_matrix_finite_nonnegative(const fc_length_qscore_2d_matrix *matrix)
{
    for (size_t i = 0; i < matrix->length_bins * matrix->qscore_bins; i++) {
        TEST_ASSERT_TRUE(isfinite(matrix->counts[i]));
        TEST_ASSERT_TRUE(matrix->counts[i] >= 0.0);
    }
}

static fc_length_qscore_2d *create_model(size_t barcode_count, size_t smoothing_points)
{
    return fc_length_qscore_2d_create(barcode_count, 64, 33, smoothing_points);
}

static void test_length_qscore_create_rejects_bad_config(void)
{
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 7, 33, 1));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 1, 1));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 33, 0));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 33, 2));
}

static void test_length_qscore_empty_model_materializes_zero_matrix(void)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_EQUAL_size_t(FC_LQ2D_LENGTH_BINS, matrix.length_bins);
    TEST_ASSERT_EQUAL_size_t(FC_LQ2D_QSCORE_BINS, matrix.qscore_bins);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, matrix.length_upper);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, FC_LQ2D_QSCORE_WIDTH, matrix.qscore_width);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, matrix_sum(&matrix));
    assert_matrix_finite_nonnegative(&matrix);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

static void test_length_qscore_single_read_lands_in_expected_qscore_bin(void)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_add_read(model, 0, 100.0, 10.25f);

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, matrix_sum(&matrix));

    size_t expected_qbin = 20;
    for (size_t q = 0; q < matrix.qscore_bins; q++) {
        double expected = (q == expected_qbin) ? 1.0 : 0.0;
        TEST_ASSERT_DOUBLE_WITHIN(1e-9, expected, column_sum(&matrix, q));
    }
    assert_matrix_finite_nonnegative(&matrix);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

static void test_length_qscore_qscore_bins_clamp_to_range(void)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_add_read(model, 0, 100.0, -2.0f);
    fc_length_qscore_2d_add_read(model, 0, 120.0, NAN);
    fc_length_qscore_2d_add_read(model, 0, 140.0, 60.0f);
    fc_length_qscore_2d_add_read(model, 0, 160.0, 999.0f);

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 4.0, matrix_sum(&matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 2.0, column_sum(&matrix, 0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 2.0, column_sum(&matrix, FC_LQ2D_QSCORE_BINS - 1));

    for (size_t q = 1; q + 1 < matrix.qscore_bins; q++) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-9, 0.0, column_sum(&matrix, q));
    }
    assert_matrix_finite_nonnegative(&matrix);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

static void test_length_qscore_barcode_state_is_isolated(void)
{
    fc_length_qscore_2d *model = create_model(3, 1);
    TEST_ASSERT_NOT_NULL(model);

    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 0));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 1));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 3));

    fc_length_qscore_2d_add_read(model, 1, 100.0, 12.0f);

    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 0));
    TEST_ASSERT_TRUE(fc_length_qscore_2d_active(model, 1));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 2));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_active(model, 3));

    fc_length_qscore_2d_matrix matrix0 = { 0 };
    fc_length_qscore_2d_matrix matrix1 = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix0));
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 1, &matrix1));
    TEST_ASSERT_EQUAL_INT(-1, fc_length_qscore_2d_materialize(model, 3, &matrix1));

    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 0.0, matrix_sum(&matrix0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, matrix_sum(&matrix1));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, column_sum(&matrix1, 24));
    assert_matrix_finite_nonnegative(&matrix0);
    assert_matrix_finite_nonnegative(&matrix1);

    fc_length_qscore_2d_matrix_destroy(&matrix0);
    fc_length_qscore_2d_matrix_destroy(&matrix1);
    fc_length_qscore_2d_destroy(model);
}

static void test_length_qscore_smoothing_preserves_observed_column_totals(void)
{
    fc_length_qscore_2d *model = create_model(1, 3);
    TEST_ASSERT_NOT_NULL(model);

    for (size_t i = 0; i < 25; i++) {
        fc_length_qscore_2d_add_read(model, 0, 100.0 + (double)i, 10.1f);
    }
    for (size_t i = 0; i < 15; i++) {
        fc_length_qscore_2d_add_read(model, 0, 300.0 + (double)i, 11.1f);
    }

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 40.0, matrix_sum(&matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 25.0, column_sum(&matrix, 20));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 15.0, column_sum(&matrix, 22));

    for (size_t q = 0; q < matrix.qscore_bins; q++) {
        if (q != 20 && q != 22) {
            TEST_ASSERT_DOUBLE_WITHIN(1e-9, 0.0, column_sum(&matrix, q));
        }
    }
    assert_matrix_finite_nonnegative(&matrix);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

static void test_length_qscore_extreme_tail_does_not_force_terminal_row_spike(void)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    for (size_t i = 0; i < 100; i++) {
        fc_length_qscore_2d_add_read(model, 0, 1000.0 + (double)i, 20.0f);
    }
    fc_length_qscore_2d_add_read(model, 0, 1000000.0, 20.0f);

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 101.0, matrix_sum(&matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 101.0, column_sum(&matrix, 40));

    double terminal = matrix.counts[40 * matrix.length_bins + matrix.length_bins - 1];
    TEST_ASSERT_TRUE(terminal < 10.0);
    assert_matrix_finite_nonnegative(&matrix);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

void length_qscore_2d_suite(void)
{
    RUN_TEST(test_length_qscore_create_rejects_bad_config);
    RUN_TEST(test_length_qscore_empty_model_materializes_zero_matrix);
    RUN_TEST(test_length_qscore_single_read_lands_in_expected_qscore_bin);
    RUN_TEST(test_length_qscore_qscore_bins_clamp_to_range);
    RUN_TEST(test_length_qscore_barcode_state_is_isolated);
    RUN_TEST(test_length_qscore_smoothing_preserves_observed_column_totals);
    RUN_TEST(test_length_qscore_extreme_tail_does_not_force_terminal_row_spike);
}
