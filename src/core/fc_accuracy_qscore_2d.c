#include "core/fc_accuracy_qscore_2d.h"

#include <math.h>
#include <stddef.h>

#include "core/fc_filter.h"
#include "core/fc_metric_qscore_2d.h"

static void aq2d_compute_range(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_range *range)
{
    (void)model;
    (void)barcode;
    range->x_lower = 0.0;
    range->x_upper = 60.0;
}

static const fc_metric_qscore_2d_spec aq2d_spec = {
    .name = "accuracy_qscore_2d",
    .x_bins = FC_AQ2D_CONCORDANCE_Q_BINS,
    .qscore_bins = FC_AQ2D_QSCORE_BINS,
    .qscore_width = FC_AQ2D_QSCORE_WIDTH,
    .client_state_size = 0,
    .compute_range = aq2d_compute_range,
};

fc_accuracy_qscore_2d *fc_accuracy_qscore_2d_create(
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
        &aq2d_spec);
}

void fc_accuracy_qscore_2d_add_read(
    fc_accuracy_qscore_2d *model,
    size_t barcode,
    const fc_read *read,
    float accuracy)
{
    if (!model || !read || !read->bam || !read->meta || !fc_is_canonical_read_record(read->bam)) {
        return;
    }
    if (!isfinite(accuracy)) {
        return;
    }
    float error_rate = (100.0f - fminf(100.0f, fmaxf(0.0f, accuracy))) / 100.0f;
    float concordance_qscore = error_rate > 0.0f
        ? -10.0f * log10f(error_rate)
        : 60.0f;
    fc_metric_qscore_2d_add_value(
        model, barcode, fminf(60.0f, concordance_qscore), read->meta->mean_q);
}
