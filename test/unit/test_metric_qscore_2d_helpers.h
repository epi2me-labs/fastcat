#ifndef FASTCAT_TEST_METRIC_QSCORE_2D_HELPERS_H
#define FASTCAT_TEST_METRIC_QSCORE_2D_HELPERS_H

#include <stddef.h>

#include "core/fc_read.h"

#define TEST_LENGTH_QSCORE_READ(LENGTH, MEAN_Q) \
    (&(fc_read) { \
        .bam = &(bam1_t) { .core = { .l_qseq = (int32_t)(LENGTH) } }, \
        .meta = &(read_meta) { .mean_q = (MEAN_Q) } \
    })

#define TEST_SPEED_QSCORE_READ(LENGTH, DURATION_S, MEAN_Q) \
    (&(fc_read) { \
        .bam = &(bam1_t) { .core = { .l_qseq = (int32_t)(LENGTH) } }, \
        .meta = &(read_meta) { .duration_s = (DURATION_S), .mean_q = (MEAN_Q) } \
    })

#define TEST_ACCURACY_QSCORE_READ(MEAN_Q) \
    (&(fc_read) { \
        .bam = &(bam1_t) { 0 }, \
        .meta = &(read_meta) { .mean_q = (MEAN_Q) } \
    })

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
