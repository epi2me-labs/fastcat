#include <math.h>

#include "unity.h"
#include "main.h"
#include "core/fc_accuracy_qscore_2d.h"
#include "test_metric_qscore_2d_helpers.h"

FC_TEST(accuracy_qscore_2d_suite, test_accuracy_qscore_bounds_qscore_bins_and_barcode_isolation)
{
    fc_accuracy_qscore_2d *model = fc_accuracy_qscore_2d_create(2, 64, 33, 3);
    TEST_ASSERT_NOT_NULL(model);

    fc_accuracy_qscore_2d_add_read(model, 0, TEST_ACCURACY_QSCORE_READ(-1.0f), -1.0f);
    fc_accuracy_qscore_2d_add_read(model, 0, TEST_ACCURACY_QSCORE_READ(60.0f), 101.0f);
    fc_accuracy_qscore_2d_add_read(model, 0, TEST_ACCURACY_QSCORE_READ(20.0f), NAN);
    fc_accuracy_qscore_2d_add_read(model, 1, TEST_ACCURACY_QSCORE_READ(10.25f), 90.0f);

    TEST_ASSERT_TRUE(fc_accuracy_qscore_2d_seen(model, 0));
    TEST_ASSERT_TRUE(fc_accuracy_qscore_2d_seen(model, 1));

    fc_accuracy_qscore_2d_matrix first = { 0 };
    fc_accuracy_qscore_2d_matrix second = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_accuracy_qscore_2d_materialize(model, 0, &first));
    TEST_ASSERT_EQUAL_INT(0, fc_accuracy_qscore_2d_materialize(model, 1, &second));
    TEST_ASSERT_EQUAL_size_t(FC_AQ2D_CONCORDANCE_Q_BINS, first.x_bins);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 60.0, first.x_upper);
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6, 2.0,
        test_metric_qscore_2d_matrix_sum(first.counts, first.x_bins, first.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6, 1.0,
        test_metric_qscore_2d_matrix_sum(second.counts, second.x_bins, second.qscore_bins));
    TEST_ASSERT_DOUBLE_WITHIN(
        1e-6, 1.0,
        test_metric_qscore_2d_column_sum(second.counts, second.x_bins, 20));
    test_metric_qscore_2d_assert_finite_nonnegative(
        first.counts, first.x_bins, first.qscore_bins);
    test_metric_qscore_2d_assert_finite_nonnegative(
        second.counts, second.x_bins, second.qscore_bins);

    fc_accuracy_qscore_2d_matrix_destroy(&first);
    fc_accuracy_qscore_2d_matrix_destroy(&second);
    fc_accuracy_qscore_2d_destroy(model);
}
