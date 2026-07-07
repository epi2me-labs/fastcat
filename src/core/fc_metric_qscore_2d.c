#include "core/fc_metric_qscore_2d.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_kll.h"
#include "core/fc_log.h"

struct fc_mq2d_barcode_state {
    fc_kll_sketch *global_x;
    fc_kll_sketch **qscore_x;
    size_t *qscore_counts;
    bool seen;
};

// generate a 2D Gaussian kernel for smoothing
static void mq2d_gaussian_kernel_2d(size_t window, double *kernel)
{
    if (window <= 1) {
        kernel[0] = 1.0;
        return;
    }

    size_t radius = window / 2;
    double sigma = (double)radius / 2.0;
    if (!(sigma > 0.0)) {
        sigma = 0.5;
    }
    double two_sigma2 = 2.0 * sigma * sigma;
    double ksum = 0.0;

    for (size_t y = 0; y < window; y++) {
        ptrdiff_t dy = (ptrdiff_t)y - (ptrdiff_t)radius;
        for (size_t x = 0; x < window; x++) {
            ptrdiff_t dx = (ptrdiff_t)x - (ptrdiff_t)radius;
            double r2 = (double)(dx * dx + dy * dy);
            double w = exp(-r2 / two_sigma2);
            kernel[y * window + x] = w;
            ksum += w;
        }
    }

    if (ksum > 0.0) {
        double inv = 1.0 / ksum;
        for (size_t i = 0; i < window * window; i++) {
            kernel[i] *= inv;
        }
    }
}

// apply direct 2D Gaussian convolution to an x by qscore matrix
static void mq2d_smooth_gaussian_2d(
    const double *in,
    size_t n_x,
    size_t n_q,
    size_t window,
    double *out)
{
    if (window <= 1 || n_x == 0 || n_q == 0) {
        memcpy(out, in, n_x * n_q * sizeof(double));
        return;
    }

    size_t radius = window / 2;
    double *kernel = fc_alloc(window * window, sizeof(double), "metric_qscore_2d gaussian kernel");
    mq2d_gaussian_kernel_2d(window, kernel);

    for (size_t q = 0; q < n_q; q++) {
        for (size_t x = 0; x < n_x; x++) {
            double sum = 0.0;
            double wsum = 0.0;

            // trim kernel footprint at matrix boundaries and track effective weight
            size_t ky_start = (q < radius) ? (radius - q) : 0;
            size_t ky_end = window - 1;
            if (q + radius >= n_q) {
                ky_end = (n_q - 1) - q + radius;
            }
            size_t kx_start = (x < radius) ? (radius - x) : 0;
            size_t kx_end = window - 1;
            if (x + radius >= n_x) {
                kx_end = (n_x - 1) - x + radius;
            }

            for (size_t ky = ky_start; ky <= ky_end; ky++) {
                size_t src_q = q + ky - radius;
                for (size_t kx = kx_start; kx <= kx_end; kx++) {
                    size_t src_x = x + kx - radius;
                    double w = kernel[ky * window + kx];
                    sum += in[src_q * n_x + src_x] * w;
                    wsum += w;
                }
            }

            out[q * n_x + x] = (wsum > 0.0) ? (sum / wsum) : in[q * n_x + x];
        }
    }

    free(kernel);
}

static void mq2d_set_empty_profile(double *rows, size_t n_x)
{
    memset(rows, 0, n_x * sizeof(double));
}

// rescale each qscore column so total mass matches observed read counts
static void mq2d_renormalize_columns_to_counts(
    double *matrix,
    size_t n_x,
    size_t n_q,
    const double *column_counts)
{
    if (!matrix || !column_counts || n_x == 0 || n_q == 0) {
        return;
    }

    for (size_t q = 0; q < n_q; q++) {
        double *col = matrix + q * n_x;
        double target = column_counts[q];

        if (target <= 0.0) {
            memset(col, 0, n_x * sizeof(double));
            continue;
        }

        double sum = 0.0;
        for (size_t i = 0; i < n_x; i++) {
            if (!isfinite(col[i]) || col[i] < 0.0) {
                col[i] = 0.0;
            }
            sum += col[i];
        }

        if (!(sum > 0.0) || !isfinite(sum)) {
            memset(col, 0, n_x * sizeof(double));
            continue;
        }

        double scale = target / sum;
        double scaled_sum = 0.0;
        for (size_t i = 0; i < n_x; i++) {
            col[i] *= scale;
            scaled_sum += col[i];
        }

        double correction = target - scaled_sum;
        if (correction != 0.0) {
            size_t best_i = 0;
            double best = col[0];
            for (size_t i = 1; i < n_x; i++) {
                if (col[i] > best) {
                    best = col[i];
                    best_i = i;
                }
            }
            col[best_i] += correction;
            if (col[best_i] < 0.0) {
                col[best_i] = 0.0;
            }
        }
    }
}

// project one constant-density interval onto fixed x rows inside the plot range
static void mq2d_add_interval_mass(
    double *rows,
    size_t n_x,
    const fc_metric_qscore_2d_range *range,
    double left,
    double right,
    double density)
{
    if (!(right > left) || !(density > 0.0) || !isfinite(density) || n_x == 0) {
        return;
    }

    if (!(range->x_upper > range->x_lower)
        || !isfinite(range->x_lower)
        || !isfinite(range->x_upper)) {
        return;
    }

    double width = range->x_upper - range->x_lower;
    double row_width = width / (double)n_x;

    if (right <= range->x_lower || left >= range->x_upper) {
        return;
    }

    if (left < range->x_lower) {
        left = range->x_lower;
    }
    if (right > range->x_upper) {
        right = range->x_upper;
    }
    if (!(right > left)) {
        return;
    }

    size_t row = (size_t)((left - range->x_lower) / row_width);
    if (row >= n_x) {
        row = n_x - 1;
    }

    while (left < right && row < n_x) {
        double row_start = range->x_lower + (double)row * row_width;
        double row_end = (row == n_x - 1) ? range->x_upper : row_start + row_width;
        double seg_end = right < row_end ? right : row_end;

        if (seg_end > left) {
            rows[row] += density * (seg_end - left);
        }

        left = seg_end;
        if (left >= right) {
            break;
        }
        row++;
    }
}

// convert irregular KLL density support into a normalized fixed-row profile
static void mq2d_project_density(
    const fc_kll_density *density,
    size_t n_x,
    const fc_metric_qscore_2d_range *range,
    double *rows)
{
    mq2d_set_empty_profile(rows, n_x);

    if (!(range->x_upper > range->x_lower)
        || !isfinite(range->x_lower)
        || !isfinite(range->x_upper)
        || !density
        || density->n == 0) {
        return;
    }

    if (density->n == 1) {
        double x = density->support[0];
        if (!isfinite(x)) {
            return;
        }
        if (x < range->x_lower || x > range->x_upper) {
            return;
        }

        double row_width = (range->x_upper - range->x_lower) / (double)n_x;
        size_t row = (size_t)((x - range->x_lower) / row_width);
        if (row >= n_x) {
            row = n_x - 1;
        }
        rows[row] = 1.0;
        return;
    }

    // reconstruct interval boundaries from support midpoints and accumulate mass
    double kept_mass = 0.0;
    for (size_t i = 0; i < density->n; i++) {
        double left = (i == 0)
            ? density->support[0] - 0.5 * (density->support[1] - density->support[0])
            : 0.5 * (density->support[i - 1] + density->support[i]);
        double right = (i + 1 < density->n)
            ? 0.5 * (density->support[i] + density->support[i + 1])
            : density->support[i] + 0.5 * (density->support[i] - density->support[i - 1]);
        mq2d_add_interval_mass(rows, n_x, range, left, right, density->values[i]);
    }

    for (size_t i = 0; i < n_x; i++) {
        if (!isfinite(rows[i]) || rows[i] < 0.0) {
            rows[i] = 0.0;
        }
        kept_mass += rows[i];
    }

    if (!(kept_mass > 0.0) || !isfinite(kept_mass)) {
        mq2d_set_empty_profile(rows, n_x);
        return;
    }

    double inv_total = 1.0 / kept_mass;
    for (size_t i = 0; i < n_x; i++) {
        rows[i] *= inv_total;
    }
}

// build one qscore column by estimating metric density and scaling to read count
static void mq2d_fill_column(
    const fc_metric_qscore_2d *model,
    const fc_mq2d_barcode_state *state,
    size_t qbin,
    const fc_metric_qscore_2d_range *range,
    double *col_counts)
{
    mq2d_set_empty_profile(col_counts, model->spec.x_bins);

    size_t observed = state->qscore_counts[qbin];
    if (observed == 0) {
        return;
    }

    fc_kll_density density = { 0 };
    int rc = -1;
    if (state->qscore_x[qbin]) {
        rc = fc_kll_compute_density_from_quantiles(
            state->qscore_x[qbin],
            model->qgrid_points,
            &density,
            true,
            model->smoothing_points);
    }
    if (rc != 0) {
        return;
    } else {
        mq2d_project_density(&density, model->spec.x_bins, range, col_counts);
        fc_kll_density_destroy(&density);
    }

    double target = (double)fc_kll_count_in_range(
        state->qscore_x[qbin],
        range->x_lower,
        range->x_upper);
    if (!(target > 0.0)) {
        mq2d_set_empty_profile(col_counts, model->spec.x_bins);
        return;
    }
    for (size_t i = 0; i < model->spec.x_bins; i++) {
        col_counts[i] *= target;
    }
}

// map mean qscore to fixed-width bins
static size_t mq2d_qscore_bin(const fc_metric_qscore_2d *model, float mean_q)
{
    if (!isfinite(mean_q)) {
        return 0;
    }

    long idx = (long)floor((double)mean_q / model->spec.qscore_width);
    if (idx < 0) {
        return 0;
    }
    if (idx >= (long)model->spec.qscore_bins) {
        return model->spec.qscore_bins - 1;
    }
    return (size_t)idx;
}

// allocate per-barcode sketch state lazily
static fc_mq2d_barcode_state *mq2d_create_barcode_state(const fc_metric_qscore_2d *model)
{
    fc_mq2d_barcode_state *state = fc_alloc(1, sizeof(*state), "metric_qscore_2d barcode state");
    state->global_x = fc_kll_create(model->kll_k);
    state->qscore_x = fc_alloc(model->spec.qscore_bins, sizeof(fc_kll_sketch *), "metric_qscore_2d qscore sketches");
    state->qscore_counts = fc_alloc(model->spec.qscore_bins, sizeof(*state->qscore_counts), "metric_qscore_2d qscore counts");
    return state;
}

// ensure barcode state exists and return it
static fc_mq2d_barcode_state *mq2d_ensure_barcode_state(fc_metric_qscore_2d *model, size_t barcode)
{
    if (barcode >= model->barcode_count) {
        FC_LOG_FATAL(
            "Internal error: barcode index %zu out of range for %s (count=%zu).",
            barcode,
            model->spec.name,
            model->barcode_count);
    }

    if (!model->barcode_states[barcode]) {
        model->barcode_states[barcode] = mq2d_create_barcode_state(model);
    }
    return model->barcode_states[barcode];
}

// free all sketches held for one barcode
static void mq2d_destroy_barcode_state(const fc_metric_qscore_2d *model, fc_mq2d_barcode_state *state)
{
    if (!state) {
        return;
    }

    fc_kll_destroy(state->global_x);
    for (size_t i = 0; i < model->spec.qscore_bins; i++) {
        fc_kll_destroy(state->qscore_x[i]);
    }
    free(state->qscore_x);
    free(state->qscore_counts);
    free(state);
}

fc_metric_qscore_2d *fc_metric_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points,
    const fc_metric_qscore_2d_spec *spec)
{
    if (!spec
        || !spec->name
        || spec->x_bins == 0
        || spec->qscore_bins == 0
        || !(spec->qscore_width > 0.0)
        || !spec->compute_range
        || kll_k < 8
        || qgrid_points < 2
        || smoothing_points == 0
        || (smoothing_points % 2) == 0) {
        return NULL;
    }

    fc_metric_qscore_2d *model = fc_alloc(1, sizeof(*model), "metric_qscore_2d model");
    model->barcode_count = barcode_count ? barcode_count : 1;
    model->kll_k = kll_k;
    model->qgrid_points = qgrid_points;
    model->smoothing_points = smoothing_points;
    model->spec = *spec;
    if (spec->client_state_size > 0) {
        model->client_state = fc_alloc(1, spec->client_state_size, "metric_qscore_2d client state");
    }
    model->barcode_states = fc_alloc(
        model->barcode_count,
        sizeof(fc_mq2d_barcode_state *),
        "metric_qscore_2d barcode state array");
    return model;
}

void fc_metric_qscore_2d_destroy(fc_metric_qscore_2d *model)
{
    if (!model) {
        return;
    }

    for (size_t i = 0; i < model->barcode_count; i++) {
        mq2d_destroy_barcode_state(model, model->barcode_states[i]);
    }
    free(model->client_state);
    free(model->barcode_states);
    free(model);
}

void fc_metric_qscore_2d_add_value(
    fc_metric_qscore_2d *model,
    size_t barcode,
    double metric_value,
    float mean_q)
{
    if (!model || !isfinite(metric_value)) {
        return;
    }

    fc_mq2d_barcode_state *state = mq2d_ensure_barcode_state(model, barcode);
    size_t qbin = mq2d_qscore_bin(model, mean_q);
    fc_kll_add(state->global_x, metric_value);

    if (!state->qscore_x[qbin]) {
        state->qscore_x[qbin] = fc_kll_create(model->kll_k);
    }
    fc_kll_add(state->qscore_x[qbin], metric_value);
    state->qscore_counts[qbin]++;
    state->seen = true;
}

void fc_metric_qscore_2d_mark_seen(fc_metric_qscore_2d *model, size_t barcode)
{
    if (!model) {
        return;
    }
    mq2d_ensure_barcode_state(model, barcode)->seen = true;
}

bool fc_metric_qscore_2d_seen(const fc_metric_qscore_2d *model, size_t barcode)
{
    if (!model || barcode >= model->barcode_count || !model->barcode_states[barcode]) {
        return false;
    }
    return model->barcode_states[barcode]->seen;
}

double fc_metric_qscore_2d_global_quantile(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    double p)
{
    if (!model || barcode >= model->barcode_count) {
        return NAN;
    }

    const fc_mq2d_barcode_state *state = model->barcode_states[barcode];
    if (!state || !state->global_x || fc_kll_count(state->global_x) == 0) {
        return NAN;
    }

    return fc_kll_compute_quantile(state->global_x, p);
}

void fc_metric_qscore_2d_compute_padded_quantile_range(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    double lower_p,
    double upper_p,
    double pad_fraction,
    double fallback_lower,
    double fallback_upper,
    fc_metric_qscore_2d_range *out)
{
    if (!out) {
        return;
    }

    out->x_lower = fallback_lower;
    out->x_upper = fallback_upper;

    if (!model
        || !(lower_p >= 0.0 && lower_p <= 1.0)
        || !(upper_p >= 0.0 && upper_p <= 1.0)
        || upper_p <= lower_p
        || !isfinite(pad_fraction)
        || pad_fraction < 0.0) {
        return;
    }

    double qlo = fc_metric_qscore_2d_global_quantile(model, barcode, lower_p);
    double qhi = fc_metric_qscore_2d_global_quantile(model, barcode, upper_p);
    if (isfinite(qlo) && isfinite(qhi) && qhi > qlo) {
        double pad = pad_fraction * (qhi - qlo);
        out->x_lower = qlo - pad;
        out->x_upper = qhi + pad;
        if (!isfinite(out->x_lower)) {
            out->x_lower = qlo;
        }
        if (!isfinite(out->x_upper)) {
            out->x_upper = qhi;
        }
        if (out->x_lower < 0.0) {
            out->x_lower = 0.0;
        }
    } else if (isfinite(qhi) && qhi > fallback_lower) {
        out->x_lower = fallback_lower;
        out->x_upper = qhi;
    }

    if (!isfinite(out->x_lower)
        || !isfinite(out->x_upper)
        || out->x_upper <= out->x_lower) {
        out->x_lower = fallback_lower;
        out->x_upper = fallback_upper;
    }
}

int fc_metric_qscore_2d_materialize(
    const fc_metric_qscore_2d *model,
    size_t barcode,
    fc_metric_qscore_2d_matrix *out)
{
    if (!model || !out || barcode >= model->barcode_count) {
        return -1;
    }

    fc_metric_qscore_2d_range range = { 0 };
    model->spec.compute_range(model, barcode, &range);

    memset(out, 0, sizeof(*out));
    out->x_bins = model->spec.x_bins;
    out->qscore_bins = model->spec.qscore_bins;
    out->x_lower = range.x_lower;
    out->x_upper = range.x_upper;
    out->qscore_width = model->spec.qscore_width;

    size_t total = model->spec.x_bins * model->spec.qscore_bins;
    out->counts = fc_alloc(total, sizeof(*out->counts), "metric_qscore_2d matrix");
    double *column_targets = fc_alloc(model->spec.qscore_bins, sizeof(*column_targets), "metric_qscore_2d column targets");

    const fc_mq2d_barcode_state *state = model->barcode_states[barcode];
    if (!state) {
        free(column_targets);
        return 0;
    }

    // build per-qscore metric profiles before optional 2D smoothing
    for (size_t q = 0; q < model->spec.qscore_bins; q++) {
        double *col = out->counts + q * model->spec.x_bins;
        mq2d_fill_column(model, state, q, &range, col);
        for (size_t i = 0; i < model->spec.x_bins; i++) {
            column_targets[q] += col[i];
        }
    }

    if (model->smoothing_points > 1) {
        size_t window = model->smoothing_points;
        if ((window % 2) == 0) {
            window += 1;
        }
        double *tmp = fc_alloc(total, sizeof(*tmp), "metric_qscore_2d smooth temp");
        mq2d_smooth_gaussian_2d(out->counts, model->spec.x_bins, model->spec.qscore_bins, window, tmp);
        memcpy(out->counts, tmp, total * sizeof(*out->counts));
        free(tmp);

        for (size_t i = 0; i < total; i++) {
            if (!isfinite(out->counts[i]) || out->counts[i] < 0.0) {
                out->counts[i] = 0.0;
            }
        }
    }

    // preserve observed qscore-bin totals after 2D smoothing
    mq2d_renormalize_columns_to_counts(
        out->counts,
        model->spec.x_bins,
        model->spec.qscore_bins,
        column_targets);
    free(column_targets);
    return 0;
}

void fc_metric_qscore_2d_matrix_destroy(fc_metric_qscore_2d_matrix *matrix)
{
    if (!matrix) {
        return;
    }
    free(matrix->counts);
    matrix->counts = NULL;
    matrix->x_bins = 0;
    matrix->qscore_bins = 0;
    matrix->x_lower = 0.0;
    matrix->x_upper = 0.0;
    matrix->qscore_width = 0.0;
}
