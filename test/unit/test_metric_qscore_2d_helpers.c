#include <math.h>
#include <stddef.h>

#include "unity.h"

#include "test_metric_qscore_2d_helpers.h"

double test_metric_qscore_2d_matrix_sum(
    const double *counts,
    size_t x_bins,
    size_t qscore_bins)
{
    double sum = 0.0;
    for (size_t i = 0; i < x_bins * qscore_bins; i++) {
        sum += counts[i];
    }
    return sum;
}

double test_metric_qscore_2d_column_sum(
    const double *counts,
    size_t x_bins,
    size_t qbin)
{
    double sum = 0.0;
    const double *col = counts + qbin * x_bins;
    for (size_t i = 0; i < x_bins; i++) {
        sum += col[i];
    }
    return sum;
}

void test_metric_qscore_2d_assert_finite_nonnegative(
    const double *counts,
    size_t x_bins,
    size_t qscore_bins)
{
    for (size_t i = 0; i < x_bins * qscore_bins; i++) {
        TEST_ASSERT_TRUE(isfinite(counts[i]));
        TEST_ASSERT_TRUE(counts[i] >= 0.0);
    }
}
