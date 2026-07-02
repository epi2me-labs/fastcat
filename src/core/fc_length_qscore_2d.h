#ifndef FASTCAT_FC_LENGTH_QSCORE_2D_H
#define FASTCAT_FC_LENGTH_QSCORE_2D_H

#include <stdbool.h>
#include <stddef.h>

#define FC_LQ2D_QSCORE_BINS 120
#define FC_LQ2D_LENGTH_BINS 1000
#define FC_LQ2D_QSCORE_WIDTH 0.5

typedef struct fc_length_qscore_2d fc_length_qscore_2d;

typedef struct {
    size_t kll_k;
    size_t qgrid_points;
    size_t smoothing_points;
} fc_length_qscore_2d_config;

typedef struct {
    size_t length_bins;
    size_t qscore_bins;
    double length_upper;
    double qscore_width;
    double *counts;
} fc_length_qscore_2d_matrix;

fc_length_qscore_2d *fc_length_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points);

void fc_length_qscore_2d_destroy(fc_length_qscore_2d *model);

void fc_length_qscore_2d_add_read(
    fc_length_qscore_2d *model,
    size_t barcode,
    double length,
    float mean_q);

bool fc_length_qscore_2d_active(const fc_length_qscore_2d *model, size_t barcode);

int fc_length_qscore_2d_materialize(
    const fc_length_qscore_2d *model,
    size_t barcode,
    fc_length_qscore_2d_matrix *out);

void fc_length_qscore_2d_matrix_destroy(fc_length_qscore_2d_matrix *matrix);

#endif
