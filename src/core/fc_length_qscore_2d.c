#include "core/fc_length_qscore_2d.h"

#include <math.h>
#include <stddef.h>

#include "core/fc_metric_qscore_2d.h"

static void lq2d_compute_range(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_range *range)
{
    range->x_lower = 0.0;
    range->x_upper = 1.0;

    double q95 = fc_metric_qscore_2d_global_quantile(model, barcode, 0.95);
    if (isfinite(q95) && q95 > 0.0) {
        range->x_upper = q95;
    }
}

static const fc_metric_qscore_2d_spec lq2d_spec = {
    .name = "length_qscore_2d",
    .x_bins = FC_LQ2D_LENGTH_BINS,
    .qscore_bins = FC_LQ2D_QSCORE_BINS,
    .qscore_width = FC_LQ2D_QSCORE_WIDTH,
    .client_state_size = 0,
    .compute_range = lq2d_compute_range,
};

fc_length_qscore_2d *fc_length_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points)
{
    return fc_metric_qscore_2d_create(
        barcode_count,
        kll_k,
        qgrid_points,
        smoothing_points,
        &lq2d_spec);
}

void fc_length_qscore_2d_add_read(
    fc_length_qscore_2d *model,
    size_t barcode,
    double length,
    float mean_q)
{
    if (!isfinite(length)) {
        return;
    }
    fc_metric_qscore_2d_add_value(model, barcode, length, mean_q);
}
