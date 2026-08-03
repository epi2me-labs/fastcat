#include <math.h>
#include <stddef.h>

#include "unity.h"
#include "main.h"
#include "core/fc_length_qscore_2d.h"
#include "test_metric_qscore_2d_helpers.h"

static fc_length_qscore_2d *create_model(size_t barcode_count, size_t smoothing_points)
{
    return fc_length_qscore_2d_create(barcode_count, 64, 33, smoothing_points);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_create_rejects_bad_config)
{
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 7, 33, 1));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 1, 1));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 33, 0));
    TEST_ASSERT_NULL(fc_length_qscore_2d_create(1, 64, 33, 2));
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_empty_model_materializes_zero_matrix)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_EQUAL_size_t(FC_LQ2D_LENGTH_BINS, matrix.x_bins);
    TEST_ASSERT_EQUAL_size_t(FC_LQ2D_QSCORE_BINS, matrix.qscore_bins);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, matrix.x_upper);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, FC_LQ2D_QSCORE_WIDTH, matrix.qscore_width);
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-12,
        0.0,
        test_metric_qscore_2d_matrix_sum(
            matrix.counts,
            matrix.x_bins,
            matrix.qscore_bins));
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_single_read_lands_in_expected_qscore_bin)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, 10.25f));

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        1.0,
        test_metric_qscore_2d_matrix_sum(
            matrix.counts,
            matrix.x_bins,
            matrix.qscore_bins));

    size_t expected_qbin = 20;
    for (size_t q = 0; q < matrix.qscore_bins; q++) {
        double expected = (q == expected_qbin) ? 1.0 : 0.0;
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-9,
            expected,
            test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, q));
    }
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_qscore_bins_clamp_to_range)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, -2.0f));
    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, NAN));
    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, 60.0f));
    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, 999.0f));

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        4.0,
        test_metric_qscore_2d_matrix_sum(
            matrix.counts,
            matrix.x_bins,
            matrix.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        2.0,
        test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, 0));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        2.0,
        test_metric_qscore_2d_column_sum(
            matrix.counts,
            matrix.x_bins,
            FC_LQ2D_QSCORE_BINS - 1));

    for (size_t q = 1; q + 1 < matrix.qscore_bins; q++) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-9,
            0.0,
            test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, q));
    }
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_barcode_state_is_isolated)
{
    fc_length_qscore_2d *model = create_model(3, 1);
    TEST_ASSERT_NOT_NULL(model);

    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 0));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 1));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 3));

    fc_length_qscore_2d_add_read(model, 1, TEST_LENGTH_QSCORE_READ(100.0, 12.0f));

    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 0));
    TEST_ASSERT_TRUE(fc_length_qscore_2d_seen(model, 1));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 2));
    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(model, 3));

    fc_length_qscore_2d_matrix matrix0 = { 0 };
    fc_length_qscore_2d_matrix matrix1 = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix0));
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 1, &matrix1));
    TEST_ASSERT_EQUAL_INT(-1, fc_length_qscore_2d_materialize(model, 3, &matrix1));

    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        0.0,
        test_metric_qscore_2d_matrix_sum(
            matrix0.counts,
            matrix0.x_bins,
            matrix0.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        1.0,
        test_metric_qscore_2d_matrix_sum(
            matrix1.counts,
            matrix1.x_bins,
            matrix1.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-9,
        1.0,
        test_metric_qscore_2d_column_sum(matrix1.counts, matrix1.x_bins, 24));
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix0.counts,
        matrix0.x_bins,
        matrix0.qscore_bins);
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix1.counts,
        matrix1.x_bins,
        matrix1.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix0);
    fc_length_qscore_2d_matrix_destroy(&matrix1);
    fc_length_qscore_2d_destroy(model);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_smoothing_preserves_observed_column_totals)
{
    fc_length_qscore_2d *model = create_model(1, 3);
    TEST_ASSERT_NOT_NULL(model);

    for (size_t i = 0; i < 25; i++) {
        fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(100.0, 10.1f));
    }
    for (size_t i = 0; i < 15; i++) {
        fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(300.0, 11.1f));
    }

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6,
        40.0,
        test_metric_qscore_2d_matrix_sum(
            matrix.counts,
            matrix.x_bins,
            matrix.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6,
        25.0,
        test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, 20));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6,
        15.0,
        test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, 22));

    for (size_t q = 0; q < matrix.qscore_bins; q++) {
        if (q != 20 && q != 22) {
            TEST_ASSERT_DOUBLE_WITHIN(
                1e-9,
                0.0,
                test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, q));
        }
    }
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}

FC_TEST(length_qscore_2d_suite, test_length_qscore_extreme_tail_does_not_force_terminal_row_spike)
{
    fc_length_qscore_2d *model = create_model(1, 1);
    TEST_ASSERT_NOT_NULL(model);

    for (size_t i = 0; i < 100; i++) {
        fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(1000.0 + (double)i, 20.0f));
    }
    fc_length_qscore_2d_add_read(model, 0, TEST_LENGTH_QSCORE_READ(1000000.0, 20.0f));

    fc_length_qscore_2d_matrix matrix = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_length_qscore_2d_materialize(model, 0, &matrix));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, matrix.x_lower);
    TEST_ASSERT_TRUE(matrix.x_upper > 1090.0);
    TEST_ASSERT_TRUE(matrix.x_upper < 1100.0);
    double total = test_metric_qscore_2d_matrix_sum(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);
    TEST_ASSERT_TRUE(total < 101.0);
    TEST_ASSERT_TRUE(total > 95.0);
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6,
        total,
        test_metric_qscore_2d_column_sum(matrix.counts, matrix.x_bins, 40));

    double terminal = matrix.counts[40 * matrix.x_bins + matrix.x_bins - 1];
    TEST_ASSERT_TRUE(terminal < 10.0);
    test_metric_qscore_2d_assert_finite_nonnegative(
        matrix.counts,
        matrix.x_bins,
        matrix.qscore_bins);

    fc_length_qscore_2d_matrix_destroy(&matrix);
    fc_length_qscore_2d_destroy(model);
}
