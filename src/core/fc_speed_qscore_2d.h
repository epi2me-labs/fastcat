#ifndef FASTCAT_FC_SPEED_QSCORE_2D_H
#define FASTCAT_FC_SPEED_QSCORE_2D_H

#include <stddef.h>

#include "core/fc_metric_qscore_2d.h"

#define FC_QTS2D_QSCORE_BINS 120
#define FC_QTS2D_SPEED_BINS 1000
#define FC_QTS2D_QSCORE_WIDTH 0.5

typedef fc_metric_qscore_2d fc_speed_qscore_2d;

typedef struct {
    size_t kll_k;
    size_t qgrid_points;
    size_t smoothing_points;
} fc_speed_qscore_2d_config;

typedef struct {
    size_t duration_missing;
    size_t duration_invalid;
} fc_speed_qscore_2d_skip_counts;

typedef fc_metric_qscore_2d_matrix fc_speed_qscore_2d_matrix;

fc_speed_qscore_2d *fc_speed_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points);

void fc_speed_qscore_2d_add_read(
    fc_speed_qscore_2d *model,
    size_t barcode,
    double length,
    double duration_s,
    float mean_q);

fc_speed_qscore_2d_skip_counts fc_speed_qscore_2d_skipped_counts(
    const fc_speed_qscore_2d *model);

#define fc_speed_qscore_2d_destroy fc_metric_qscore_2d_destroy
#define fc_speed_qscore_2d_seen fc_metric_qscore_2d_seen
#define fc_speed_qscore_2d_materialize fc_metric_qscore_2d_materialize
#define fc_speed_qscore_2d_matrix_destroy fc_metric_qscore_2d_matrix_destroy

#endif
