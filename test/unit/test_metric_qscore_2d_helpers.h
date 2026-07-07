#ifndef FASTCAT_TEST_METRIC_QSCORE_2D_HELPERS_H
#define FASTCAT_TEST_METRIC_QSCORE_2D_HELPERS_H

#include <stddef.h>

double test_metric_qscore_2d_matrix_sum(
    const double *counts,
    size_t x_bins,
    size_t qscore_bins);

double test_metric_qscore_2d_column_sum(
    const double *counts,
    size_t x_bins,
    size_t qbin);

void test_metric_qscore_2d_assert_finite_nonnegative(
    const double *counts,
    size_t x_bins,
    size_t qscore_bins);

#endif
