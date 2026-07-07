#ifndef FASTCAT_FC_METRIC_QSCORE_2D_H
#define FASTCAT_FC_METRIC_QSCORE_2D_H

#include <stdbool.h>
#include <stddef.h>

typedef struct fc_mq2d_barcode_state fc_mq2d_barcode_state;
typedef struct fc_metric_qscore_2d fc_metric_qscore_2d;

typedef struct {
    double x_lower;
    double x_upper;
} fc_metric_qscore_2d_range;

typedef void (*fc_metric_qscore_2d_range_fn)(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_range *out);

typedef struct {
    const char *name;
    size_t x_bins;
    size_t qscore_bins;
    double qscore_width;
    size_t client_state_size;
    fc_metric_qscore_2d_range_fn compute_range;
} fc_metric_qscore_2d_spec;

typedef struct {
    size_t x_bins;
    size_t qscore_bins;
    double x_lower;
    double x_upper;
    double qscore_width;
    double *counts;
} fc_metric_qscore_2d_matrix;

struct fc_metric_qscore_2d {
    size_t barcode_count;
    size_t kll_k;
    size_t qgrid_points;
    size_t smoothing_points;
    fc_metric_qscore_2d_spec spec;
    void *client_state;
    fc_mq2d_barcode_state **barcode_states;
};

fc_metric_qscore_2d *fc_metric_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points,
    const fc_metric_qscore_2d_spec *spec);

void fc_metric_qscore_2d_destroy(fc_metric_qscore_2d *model);

void fc_metric_qscore_2d_add_value(
    fc_metric_qscore_2d *model,
    size_t barcode,
    double metric_value,
    float mean_q);

void fc_metric_qscore_2d_mark_seen(fc_metric_qscore_2d *model, size_t barcode);

bool fc_metric_qscore_2d_seen(const fc_metric_qscore_2d *model, size_t barcode);

double fc_metric_qscore_2d_global_quantile(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    double p);

void fc_metric_qscore_2d_compute_padded_quantile_range(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    double lower_p,
    double upper_p,
    double pad_fraction,
    double fallback_lower,
    double fallback_upper,
    fc_metric_qscore_2d_range *out);

int fc_metric_qscore_2d_materialize(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_matrix *out);

void fc_metric_qscore_2d_matrix_destroy(fc_metric_qscore_2d_matrix *matrix);

#endif
