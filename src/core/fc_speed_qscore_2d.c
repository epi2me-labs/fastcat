#include "core/fc_speed_qscore_2d.h"

#include <math.h>
#include <stddef.h>

#include "core/fc_metric_qscore_2d.h"

typedef struct {
    size_t skipped_duration_missing;
    size_t skipped_duration_invalid;
} qts2d_client_state;

static void qts2d_compute_range(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_range *range)
{
    fc_metric_qscore_2d_compute_padded_quantile_range(
        model,
        barcode,
        0.05,
        0.95,
        0.1,
        0.0,
        1.0,
        range);
}

static const fc_metric_qscore_2d_spec qts2d_spec = {
    .name = "speed_qscore_2d",
    .x_bins = FC_QTS2D_SPEED_BINS,
    .qscore_bins = FC_QTS2D_QSCORE_BINS,
    .qscore_width = FC_QTS2D_QSCORE_WIDTH,
    .client_state_size = sizeof(qts2d_client_state),
    .compute_range = qts2d_compute_range,
};

fc_speed_qscore_2d *fc_speed_qscore_2d_create(
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
        &qts2d_spec);
}

void fc_speed_qscore_2d_add_read(
    fc_speed_qscore_2d *model,
    size_t barcode,
    double length,
    double duration_s,
    float mean_q)
{
    fc_metric_qscore_2d_mark_seen(model, barcode);

    qts2d_client_state *state = model->client_state;
    if (!isfinite(duration_s)) {
        state->skipped_duration_missing++;
        return;
    }
    if (duration_s <= 0.0 || !isfinite(length)) {
        state->skipped_duration_invalid++;
        return;
    }

    double speed = length / duration_s;
    if (!isfinite(speed) || speed <= 0.0) {
        state->skipped_duration_invalid++;
        return;
    }

    fc_metric_qscore_2d_add_value(model, barcode, speed, mean_q);
}

fc_speed_qscore_2d_skip_counts fc_speed_qscore_2d_skipped_counts(
    const fc_speed_qscore_2d *model)
{
    fc_speed_qscore_2d_skip_counts counts = { 0 };
    const qts2d_client_state *state = model ? model->client_state : NULL;
    if (!state) {
        return counts;
    }
    counts.duration_missing = state->skipped_duration_missing;
    counts.duration_invalid = state->skipped_duration_invalid;
    return counts;
}
