#include "core/fc_length_qscore_2d.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_kll.h"
#include "core/fc_log.h"

typedef struct {
    fc_kll_sketch *global_length;
    fc_kll_sketch *qscore_lengths[FC_LQ2D_QSCORE_BINS];
    size_t qscore_counts[FC_LQ2D_QSCORE_BINS];
    bool seen;
} fc_lq2d_barcode_state;

struct fc_length_qscore_2d {
    size_t barcode_count;
    size_t kll_k;
    size_t qgrid_points;
    size_t smoothing_points;
    fc_lq2d_barcode_state **barcode_states;
};

// generate a 2D Gaussian kernel for smoothing
static void lq2d_gaussian_kernel_2d(size_t window, double *kernel)
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

// apply direct 2D gaussian convolution to a length x qscore matrix
static void lq2d_smooth_gaussian_2d(
    const double *in,
    size_t n_len,
    size_t n_q,
    size_t window,
    double *out)
{
    if (window <= 1 || n_len == 0 || n_q == 0) {
        memcpy(out, in, n_len * n_q * sizeof(double));
        return;
    }

    size_t radius = window / 2;
    double *kernel = fc_alloc(window * window, sizeof(double), "length_qscore_2d gaussian kernel");
    lq2d_gaussian_kernel_2d(window, kernel);

    for (size_t q = 0; q < n_q; q++) {
        for (size_t len = 0; len < n_len; len++) {
            double sum = 0.0;
            double wsum = 0.0;

            // trim kernel footprint at matrix boundaries and track effective weight
            size_t ky_start = (q < radius) ? (radius - q) : 0;
            size_t ky_end = window - 1;
            if (q + radius >= n_q) {
                ky_end = (n_q - 1) - q + radius;
            }
            size_t kx_start = (len < radius) ? (radius - len) : 0;
            size_t kx_end = window - 1;
            if (len + radius >= n_len) {
                kx_end = (n_len - 1) - len + radius;
            }

            for (size_t ky = ky_start; ky <= ky_end; ky++) {
                size_t src_q = q + ky - radius;
                for (size_t kx = kx_start; kx <= kx_end; kx++) {
                    size_t src_len = len + kx - radius;
                    double w = kernel[ky * window + kx];
                    sum += in[src_q * n_len + src_len] * w;
                    wsum += w;
                }
            }

            out[q * n_len + len] = (wsum > 0.0) ? (sum / wsum) : in[q * n_len + len];
        }
    }

    free(kernel);
}

// rescale each qscore column so total mass matches observed read counts
static void lq2d_renormalize_columns_to_counts(
    double *matrix,
    size_t n_len,
    size_t n_q,
    const size_t *column_counts)
{
    if (!matrix || !column_counts || n_len == 0 || n_q == 0) {
        return;
    }

    for (size_t q = 0; q < n_q; q++) {
        double *col = matrix + q * n_len;
        double target = (double)column_counts[q];

        if (target <= 0.0) {
            memset(col, 0, n_len * sizeof(double));
            continue;
        }

        double sum = 0.0;
        for (size_t i = 0; i < n_len; i++) {
            if (!isfinite(col[i]) || col[i] < 0.0) {
                col[i] = 0.0;
            }
            sum += col[i];
        }

        if (!(sum > 0.0) || !isfinite(sum)) {
            memset(col, 0, n_len * sizeof(double));
            col[n_len - 1] = target;
            continue;
        }

        double scale = target / sum;
        double scaled_sum = 0.0;
        for (size_t i = 0; i < n_len; i++) {
            col[i] *= scale;
            scaled_sum += col[i];
        }
        col[n_len - 1] += target - scaled_sum;
        if (col[n_len - 1] < 0.0) {
            col[n_len - 1] = 0.0;
        }
    }
}

// project one constant-density interval onto fixed length rows inside the plot range
static void lq2d_add_interval_mass(double *rows, double upper, double left, double right, double density)
{
    if (!(right > left) || !(density > 0.0) || !isfinite(density)) {
        return;
    }

    double row_width = upper / (double)FC_LQ2D_LENGTH_BINS;

    if (right <= 0.0) {
        return;
    }
    if (left < 0.0) {
        left = 0.0;
    }
    if (left >= upper) {
        return;
    }
    if (right > upper) {
        right = upper;
    }
    if (!(right > left)) {
        return;
    }

    size_t row = (size_t)(left / row_width);
    if (row >= FC_LQ2D_LENGTH_BINS) {
        row = FC_LQ2D_LENGTH_BINS - 1;
    }

    while (left < right && row < FC_LQ2D_LENGTH_BINS) {
        double row_start = (double)row * row_width;
        double row_end = (row == FC_LQ2D_LENGTH_BINS - 1) ? upper : row_start + row_width;
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

// convert irregular KLL density support into a normalized fixed-row length profile
static void lq2d_project_density(const fc_kll_density *density, double upper, double *rows)
{
    memset(rows, 0, FC_LQ2D_LENGTH_BINS * sizeof(double));

    if (!(upper > 0.0) || !isfinite(upper) || !density || density->n == 0) {
        rows[FC_LQ2D_LENGTH_BINS - 1] = 1.0;
        return;
    }

    if (density->n == 1) {
        double x = density->support[0];
        size_t row = FC_LQ2D_LENGTH_BINS - 1;
        if (isfinite(x) && x <= 0.0) {
            row = 0;
        } else if (isfinite(x) && x < upper) {
            double row_width = upper / (double)FC_LQ2D_LENGTH_BINS;
            row = (size_t)(x / row_width);
            if (row >= FC_LQ2D_LENGTH_BINS) {
                row = FC_LQ2D_LENGTH_BINS - 1;
            }
        }
        rows[row] = 1.0;
        return;
    }

    // reconstruct interval boundaries from support midpoints and accumulate mass
    for (size_t i = 0; i < density->n; i++) {
        double left = (i == 0)
            ? density->support[0] - 0.5 * (density->support[1] - density->support[0])
            : 0.5 * (density->support[i - 1] + density->support[i]);
        double right = (i + 1 < density->n)
            ? 0.5 * (density->support[i] + density->support[i + 1])
            : density->support[i] + 0.5 * (density->support[i] - density->support[i - 1]);
        lq2d_add_interval_mass(rows, upper, left, right, density->values[i]);
    }

    // normalize projected mass to unit area and guard numeric failures
    double total = 0.0;
    for (size_t i = 0; i < FC_LQ2D_LENGTH_BINS; i++) {
        if (!isfinite(rows[i]) || rows[i] < 0.0) {
            rows[i] = 0.0;
        }
        total += rows[i];
    }

    if (!(total > 0.0) || !isfinite(total)) {
        memset(rows, 0, FC_LQ2D_LENGTH_BINS * sizeof(double));
        rows[FC_LQ2D_LENGTH_BINS - 1] = 1.0;
        return;
    }

    double inv_total = 1.0 / total;
    for (size_t i = 0; i < FC_LQ2D_LENGTH_BINS; i++) {
        rows[i] *= inv_total;
    }
}

// build one qscore column by estimating length density and scaling to read count
static void lq2d_fill_column(
    const fc_length_qscore_2d *model,
    const fc_lq2d_barcode_state *bc,
    size_t qbin,
    double length_upper,
    double *col_counts)
{
    memset(col_counts, 0, FC_LQ2D_LENGTH_BINS * sizeof(double));

    size_t observed = bc->qscore_counts[qbin];
    if (observed == 0) {
        return;
    }

    fc_kll_density density = { 0 };
    int rc = -1;
    if (bc->qscore_lengths[qbin]) {
        rc = fc_kll_compute_density_from_quantiles(
            bc->qscore_lengths[qbin],
            model->qgrid_points,
            &density,
            true,
            model->smoothing_points);
    }
    // if density estimation fails, place all count mass in the terminal row
    if (rc != 0) {
        col_counts[FC_LQ2D_LENGTH_BINS - 1] = (double)observed;
        return;
    }

    lq2d_project_density(&density, length_upper, col_counts);
    fc_kll_density_destroy(&density);

    double col_sum = 0.0;
    for (size_t i = 0; i < FC_LQ2D_LENGTH_BINS; i++) {
        col_sum += col_counts[i];
    }

    // enforce exact column sum equal to observed qscore-bin reads
    double target = (double)observed;
    if (!(col_sum > 0.0) || !isfinite(col_sum)) {
        memset(col_counts, 0, FC_LQ2D_LENGTH_BINS * sizeof(double));
        col_counts[FC_LQ2D_LENGTH_BINS - 1] = target;
        return;
    }

    double scale = target / col_sum;
    double scaled_sum = 0.0;
    for (size_t i = 0; i < FC_LQ2D_LENGTH_BINS; i++) {
        col_counts[i] *= scale;
        scaled_sum += col_counts[i];
    }

    col_counts[FC_LQ2D_LENGTH_BINS - 1] += target - scaled_sum;
    if (col_counts[FC_LQ2D_LENGTH_BINS - 1] < 0.0) {
        col_counts[FC_LQ2D_LENGTH_BINS - 1] = 0.0;
    }
}

// map mean qscore to fixed 0.5-wide bins across [0, 60)
static size_t lq2d_qscore_bin(float mean_q)
{
    if (!isfinite(mean_q)) {
        return 0;
    }

    long idx = (long)floor((double)mean_q * 2.0);
    if (idx < 0) {
        return 0;
    }
    if (idx >= (long)FC_LQ2D_QSCORE_BINS) {
        return FC_LQ2D_QSCORE_BINS - 1;
    }
    return (size_t)idx;
}

// allocate per-barcode sketch state lazily
static fc_lq2d_barcode_state *lq2d_create_barcode_state(size_t kll_k)
{
    fc_lq2d_barcode_state *state = fc_alloc(1, sizeof(*state), "length_qscore_2d barcode state");
    state->global_length = fc_kll_create(kll_k);
    return state;
}

// ensure barcode state exists and return it
static fc_lq2d_barcode_state *lq2d_ensure_barcode_state(fc_length_qscore_2d *model, size_t barcode)
{
    if (barcode >= model->barcode_count) {
        FC_LOG_FATAL("Internal error: barcode index %zu out of range for length_qscore_2d (count=%zu).",
            barcode, model->barcode_count);
    }

    if (!model->barcode_states[barcode]) {
        model->barcode_states[barcode] = lq2d_create_barcode_state(model->kll_k);
    }
    return model->barcode_states[barcode];
}

// free all sketches held for one barcode
static void lq2d_destroy_barcode_state(fc_lq2d_barcode_state *state)
{
    if (!state) {
        return;
    }

    fc_kll_destroy(state->global_length);
    for (size_t i = 0; i < FC_LQ2D_QSCORE_BINS; i++) {
        fc_kll_destroy(state->qscore_lengths[i]);
    }
    free(state);
}

fc_length_qscore_2d *fc_length_qscore_2d_create(
    size_t barcode_count,
    size_t kll_k,
    size_t qgrid_points,
    size_t smoothing_points)
{
    if (kll_k < 8 || qgrid_points < 2 || smoothing_points == 0 || (smoothing_points % 2) == 0) {
        return NULL;
    }

    fc_length_qscore_2d *model = fc_alloc(1, sizeof(*model), "length_qscore_2d model");
    model->barcode_count = barcode_count ? barcode_count : 1;
    model->kll_k = kll_k;
    model->qgrid_points = qgrid_points;
    model->smoothing_points = smoothing_points;
    model->barcode_states = fc_alloc(
        model->barcode_count,
        sizeof(fc_lq2d_barcode_state *),
        "length_qscore_2d barcode state array");
    return model;
}

void fc_length_qscore_2d_destroy(fc_length_qscore_2d *model)
{
    if (!model) {
        return;
    }

    for (size_t i = 0; i < model->barcode_count; i++) {
        lq2d_destroy_barcode_state(model->barcode_states[i]);
    }
    free(model->barcode_states);
    free(model);
}

void fc_length_qscore_2d_add_read(
    fc_length_qscore_2d *model,
    size_t barcode,
    double length,
    float mean_q)
{
    if (!model) {
        return;
    }

    fc_lq2d_barcode_state *bc = lq2d_ensure_barcode_state(model, barcode);
    size_t qbin = lq2d_qscore_bin(mean_q);
    fc_kll_add(bc->global_length, length);

    // allocate qscore-specific sketch on first use
    if (!bc->qscore_lengths[qbin]) {
        bc->qscore_lengths[qbin] = fc_kll_create(model->kll_k);
    }
    fc_kll_add(bc->qscore_lengths[qbin], length);
    bc->qscore_counts[qbin]++;
    bc->seen = true;
}

bool fc_length_qscore_2d_seen(const fc_length_qscore_2d *model, size_t barcode)
{
    if (!model || barcode >= model->barcode_count || !model->barcode_states[barcode]) {
        return false;
    }
    return model->barcode_states[barcode]->seen;
}

int fc_length_qscore_2d_materialize(
    const fc_length_qscore_2d *model,
    size_t barcode,
    fc_length_qscore_2d_matrix *out)
{
    if (!model || !out || barcode >= model->barcode_count) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    out->length_bins = FC_LQ2D_LENGTH_BINS;
    out->qscore_bins = FC_LQ2D_QSCORE_BINS;
    out->qscore_width = FC_LQ2D_QSCORE_WIDTH;

    const fc_lq2d_barcode_state *bc = model->barcode_states[barcode];
    double length_upper = 1.0;
    if (bc && bc->global_length && fc_kll_count(bc->global_length) > 0) {
        length_upper = fc_kll_compute_quantile(bc->global_length, 0.95);
        if (!isfinite(length_upper) || length_upper <= 0.0) {
            length_upper = 1.0;
        }
    }
    out->length_upper = length_upper;

    size_t total = FC_LQ2D_LENGTH_BINS * FC_LQ2D_QSCORE_BINS;
    out->counts = fc_alloc(total, sizeof(*out->counts), "length_qscore_2d matrix");

    if (!bc) {
        return 0;
    }

    // build per-qscore length profiles before optional 2D smoothing
    for (size_t q = 0; q < FC_LQ2D_QSCORE_BINS; q++) {
        lq2d_fill_column(model, bc, q, length_upper, out->counts + q * FC_LQ2D_LENGTH_BINS);
    }

    if (model->smoothing_points > 1) {
        size_t window = model->smoothing_points;
        if ((window % 2) == 0) {
            window += 1;
        }
        double *tmp = fc_alloc(total, sizeof(*tmp), "length_qscore_2d smooth temp");
        lq2d_smooth_gaussian_2d(out->counts, FC_LQ2D_LENGTH_BINS, FC_LQ2D_QSCORE_BINS, window, tmp);
        memcpy(out->counts, tmp, total * sizeof(*out->counts));
        free(tmp);

        for (size_t i = 0; i < total; i++) {
            if (!isfinite(out->counts[i]) || out->counts[i] < 0.0) {
                out->counts[i] = 0.0;
            }
        }
    }

    // preserve observed qscore-bin totals after 2D smoothing
    lq2d_renormalize_columns_to_counts(
        out->counts,
        FC_LQ2D_LENGTH_BINS,
        FC_LQ2D_QSCORE_BINS,
        bc->qscore_counts);
    return 0;
}

void fc_length_qscore_2d_matrix_destroy(fc_length_qscore_2d_matrix *matrix)
{
    if (!matrix) {
        return;
    }
    free(matrix->counts);
    matrix->counts = NULL;
    matrix->length_bins = 0;
    matrix->qscore_bins = 0;
    matrix->length_upper = 0.0;
    matrix->qscore_width = 0.0;
}
