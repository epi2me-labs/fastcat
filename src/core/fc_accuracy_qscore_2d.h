#ifndef FASTCAT_FC_ACCURACY_QSCORE_2D_H
#define FASTCAT_FC_ACCURACY_QSCORE_2D_H

#include <stddef.h>

#include "core/fc_metric_qscore_2d.h"
#include "core/fc_read.h"

#define FC_AQ2D_QSCORE_BINS 120
#define FC_AQ2D_CONCORDANCE_Q_BINS 120
#define FC_AQ2D_QSCORE_WIDTH 0.5

typedef fc_metric_qscore_2d fc_accuracy_qscore_2d;

typedef struct {
    size_t kll_k;
    size_t qgrid_points;
    size_t smoothing_points;
} fc_accuracy_qscore_2d_config;

typedef fc_metric_qscore_2d_matrix fc_accuracy_qscore_2d_matrix;

fc_accuracy_qscore_2d *fc_accuracy_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points);

void fc_accuracy_qscore_2d_add_read(
    fc_accuracy_qscore_2d *model,
    size_t barcode,
    const fc_read *read,
    float accuracy);

#define fc_accuracy_qscore_2d_destroy fc_metric_qscore_2d_destroy
#define fc_accuracy_qscore_2d_seen fc_metric_qscore_2d_seen
#define fc_accuracy_qscore_2d_materialize fc_metric_qscore_2d_materialize
#define fc_accuracy_qscore_2d_matrix_destroy fc_metric_qscore_2d_matrix_destroy

#endif
