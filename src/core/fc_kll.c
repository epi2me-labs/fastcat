#include "core/fc_kll.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"

// KLL sketch implementation
//
// The sketch maintains a hierarchy of "compactors" (levels) where each level
// stores up to k items. When a level fills, half of its items are promoted
// to the next level (with increased weight). This creates a compressed
// representation of the data distribution.

typedef struct {
    double *items; // sorted items in this level
    size_t count; // current number of items
    size_t capacity; // k items per level
} kll_compactor;

struct fc_kll_sketch {
    size_t k; // accuracy parameter
    size_t n; // total items added
    kll_compactor *levels; // array of compactors
    size_t num_levels; // current number of levels
    size_t max_levels; // allocated capacity for levels
    uint32_t rng_state; // PRNG state for unbiased sampling during compaction
};

// temporary struct for sorting items with their weights
typedef struct {
    double item;
    size_t weight;
} item_weight_pair;

static size_t get_sorted_items(const fc_kll_sketch *sketch, double **items_out, size_t **weights_out);

// compare function for sorting item-weight pairs
static int compare_item_weight_pairs(const void *a, const void *b)
{
    const item_weight_pair *pa = (const item_weight_pair *)a;
    const item_weight_pair *pb = (const item_weight_pair *)b;
    if (pa->item < pb->item) {
        return -1;
    }
    if (pa->item > pb->item) {
        return 1;
    }
    return 0;
}

// Simple xorshift32 PRNG for generating random offset during compaction
// This ensures unbiased sampling when compacting levels
static inline uint32_t xorshift32(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

fc_kll_sketch *fc_kll_create(size_t k)
{
    if (k < 8) {
        FC_LOG_FATAL("KLL sketch k parameter must be >= 8 (got %zu)", k);
    }

    fc_kll_sketch *sketch = fc_alloc(1, sizeof(fc_kll_sketch), "kll_sketch");
    sketch->k = k;
    sketch->n = 0;
    sketch->rng_state = 0xDEADBEEF;

    // start with a reasonable number of levels
    sketch->max_levels = 16;
    sketch->num_levels = 1;
    sketch->levels = fc_alloc(sketch->max_levels, sizeof(kll_compactor), "kll_levels");

    // initialize first level
    sketch->levels[0].capacity = k;
    sketch->levels[0].count = 0;
    sketch->levels[0].items = fc_alloc(k, sizeof(double), "kll_items");

    return sketch;
}

void fc_kll_destroy(fc_kll_sketch *sketch)
{
    if (sketch == NULL) {
        return;
    }
    for (size_t i = 0; i < sketch->num_levels; i++) {
        free(sketch->levels[i].items);
    }
    free(sketch->levels);
    free(sketch);
}

// ensure we have enough levels, growing if needed
static void ensure_level_exists(fc_kll_sketch *sketch, size_t level)
{
    if (level >= sketch->max_levels) {
        size_t new_max = sketch->max_levels * 2;
        sketch->levels = fc_realloc(sketch->levels, new_max * sizeof(kll_compactor), "kll_levels");
        sketch->max_levels = new_max;
    }

    while (sketch->num_levels <= level) {
        kll_compactor *comp = &sketch->levels[sketch->num_levels];
        comp->capacity = sketch->k;
        comp->count = 0;
        comp->items = fc_alloc(sketch->k, sizeof(double), "kll_items");
        sketch->num_levels++;
    }
}

// compact a level following KLL algorithm:
// 1. sort items
// 2. if odd count, save one random endpoint to keep at this level
// 3. promote every other item (starting from random offset) to next level
// 4. clear this level and restore the saved item if any
static void compact_level(fc_kll_sketch *sketch, size_t level)
{
    kll_compactor *comp = &sketch->levels[level];

    // sort items in this level
    // if odd number of items, keep one random endpoint at this level
    qsort(comp->items, comp->count, sizeof(double), cmp_double);
    bool is_odd = (comp->count % 2) == 1;
    double saved_item = 0.0;
    size_t start_idx = 0;
    size_t count_to_process = comp->count;

    if (is_odd) {
        bool keep_first = xorshift32(&sketch->rng_state) >> 31; // random boolean
        saved_item = keep_first ? comp->items[0] : comp->items[comp->count - 1];
        start_idx = keep_first ? 1 : 0;
        count_to_process--;
    }

    ensure_level_exists(sketch, level + 1);

    // refresh pointers after potential realloc in ensure_level_exists
    comp = &sketch->levels[level];
    kll_compactor *next = &sketch->levels[level + 1];

    // promote every other item to next level
    size_t offset = (size_t)(xorshift32(&sketch->rng_state) >> 31);

    for (size_t i = start_idx + offset; i < start_idx + count_to_process; i += 2) {
        if (sketch->levels[level + 1].count >= sketch->levels[level + 1].capacity) {
            compact_level(sketch, level + 1);
            // refresh pointers after recursive compact which may realloc
            comp = &sketch->levels[level];
            next = &sketch->levels[level + 1];
        }
        next->items[next->count++] = comp->items[i];
    }

    // clear this level and restore saved item if any
    comp->count = 0;
    if (is_odd) {
        comp->items[0] = saved_item;
        comp->count = 1;
    }
}

void fc_kll_add(fc_kll_sketch *sketch, double value)
{
    if (sketch == NULL || !isfinite(value)) {
        return;
    }

    sketch->n++;
    if (sketch->levels[0].count >= sketch->levels[0].capacity) {
        compact_level(sketch, 0);
    }
    // refresh pointer after potential realloc in compact_level
    sketch->levels[0].items[sketch->levels[0].count++] = value;
}

static void add_retained_item_to_level(fc_kll_sketch *sketch, size_t level, double value)
{
    ensure_level_exists(sketch, level);

    while (sketch->levels[level].count >= sketch->levels[level].capacity) {
        compact_level(sketch, level);
        ensure_level_exists(sketch, level);
    }
    sketch->levels[level].items[sketch->levels[level].count++] = value;
}

void fc_kll_merge(fc_kll_sketch *dst, const fc_kll_sketch *src)
{
    if (dst == NULL || src == NULL || src->n == 0) {
        return;
    }
    if (dst == src) {
        FC_LOG_FATAL("Cannot merge a KLL sketch into itself.");
    }
    if (dst->k != src->k) {
        FC_LOG_FATAL("Cannot merge KLL sketches with different k values (%zu != %zu).", dst->k, src->k);
    }
    if (dst->n > SIZE_MAX - src->n) {
        FC_LOG_FATAL("KLL sketch count overflow during merge.");
    }

    for (size_t level = 0; level < src->num_levels; level++) {
        const kll_compactor *comp = &src->levels[level];
        for (size_t i = 0; i < comp->count; i++) {
            add_retained_item_to_level(dst, level, comp->items[i]);
        }
    }
    dst->n += src->n;
}

size_t fc_kll_count(const fc_kll_sketch *sketch)
{
    return sketch == NULL ? 0 : sketch->n;
}

size_t fc_kll_count_in_range(const fc_kll_sketch *sketch, double lower, double upper)
{
    if (sketch == NULL || !isfinite(lower) || !isfinite(upper) || upper < lower) {
        return 0;
    }

    double *items = NULL;
    size_t *weights = NULL;
    size_t count = get_sorted_items(sketch, &items, &weights);
    if (count == 0) {
        return 0;
    }

    size_t total = 0;
    for (size_t i = 0; i < count; i++) {
        if (items[i] >= lower && items[i] <= upper) {
            total += weights[i];
        }
    }

    free(items);
    free(weights);
    return total;
}

size_t fc_kll_memory_bytes(const fc_kll_sketch *sketch)
{
    if (sketch == NULL) {
        return 0;
    }

    size_t bytes = sizeof(fc_kll_sketch);
    bytes += sketch->max_levels * sizeof(kll_compactor);
    for (size_t i = 0; i < sketch->num_levels; i++) {
        bytes += sketch->levels[i].capacity * sizeof(double);
    }

    return bytes;
}

// merge all levels into a sorted array with weights
// returns number of unique items
static size_t get_sorted_items(const fc_kll_sketch *sketch, double **items_out, size_t **weights_out)
{
    if (sketch->n == 0) {
        *items_out = NULL;
        *weights_out = NULL;
        return 0;
    }

    // count total items across all levels
    size_t total = 0;
    for (size_t i = 0; i < sketch->num_levels; i++) {
        total += sketch->levels[i].count;
    }

    if (total == 0) {
        *items_out = NULL;
        *weights_out = NULL;
        return 0;
    }

    // collect items with their weights (2^level) into temporary array
    item_weight_pair *pairs = fc_alloc(total, sizeof(item_weight_pair), "kll_pairs_temp");

    size_t idx = 0;
    for (size_t level = 0; level < sketch->num_levels; level++) {
        size_t weight = 1ULL << level; // 2^level
        kll_compactor *comp = &sketch->levels[level];
        for (size_t i = 0; i < comp->count; i++) {
            pairs[idx].item = comp->items[i];
            pairs[idx].weight = weight;
            idx++;
        }
    }

    // sort by item value, then extract.
    qsort(pairs, total, sizeof(item_weight_pair), compare_item_weight_pairs);
    double *items = fc_alloc(total, sizeof(double), "kll_items_temp");
    size_t *weights = fc_alloc(total, sizeof(size_t), "kll_weights_temp");

    for (size_t i = 0; i < total; i++) {
        items[i] = pairs[i].item;
        weights[i] = pairs[i].weight;
    }

    free(pairs);

    *items_out = items;
    *weights_out = weights;
    return total;
}

double fc_kll_compute_quantile(const fc_kll_sketch *sketch, double p)
{
    double result = NAN;
    fc_kll_compute_quantiles(sketch, &p, 1, &result);
    return result;
}

fc_kll_quantile_septuple fc_kll_compute_quantile_septuple(const fc_kll_sketch *sketch)
{
    double p[7] = { 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95 };
    double out[7] = { NAN, NAN, NAN, NAN, NAN, NAN, NAN };
    fc_kll_quantile_septuple result;

    fc_kll_compute_quantiles(sketch, p, 7, out);

    result.q05 = out[0];
    result.q10 = out[1];
    result.q25 = out[2];
    result.q50 = out[3];
    result.q75 = out[4];
    result.q90 = out[5];
    result.q95 = out[6];

    return result;
}

// Find first index i such that cumsum[i] > order_idx (0-based order statistic)
static size_t bisect(const size_t *cumsum, size_t n, size_t order_idx)
{
    size_t left = 0;
    size_t right = n;

    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (cumsum[mid] > order_idx) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

int fc_kll_compute_quantiles(const fc_kll_sketch *sketch, const double *ps, size_t n, double *out)
{
    if (sketch == NULL || ps == NULL || out == NULL || n == 0) {
        return -1;
    }

    if (sketch->n == 0) {
        for (size_t i = 0; i < n; i++) {
            out[i] = NAN;
        }
        return 0;
    }

    double *items = NULL;
    size_t *weights = NULL;
    size_t count = get_sorted_items(sketch, &items, &weights);

    if (count == 0) {
        for (size_t i = 0; i < n; i++) {
            out[i] = NAN;
        }
        return 0;
    }

    // build cumulative weight array whilst getting total
    size_t total_weight = 0;
    size_t *cumsum = fc_alloc(count, sizeof(size_t), "kll_cumsum");
    for (size_t i = 0; i < count; i++) {
        total_weight += weights[i];
        cumsum[i] = total_weight;
    }

    if (total_weight == 0) {
        for (size_t i = 0; i < n; i++) {
            out[i] = NAN;
        }
        free(cumsum);
        free(items);
        free(weights);
        return 0;
    }

    // interpolate quantiles between neighbouring items taking into
    // account their weights.
    // t = p * (W - 1); f = fract(t)
    // x = (1-f) * x[floor(t)] + f * x[ceil(t)]
    for (size_t q = 0; q < n; q++) {
        if (ps[q] < 0.0 || ps[q] > 1.0) {
            out[q] = NAN;
            continue;
        }

        if (total_weight == 1 || ps[q] <= 0.0) {
            out[q] = items[0];
            continue;
        }

        if (ps[q] >= 1.0) {
            out[q] = items[count - 1];
            continue;
        }

        double t = ps[q] * (double)(total_weight - 1);
        size_t idx_lo = (size_t)t;
        size_t idx_hi = idx_lo + (t > (double)idx_lo ? 1 : 0);
        double frac = t - (double)idx_lo;

        size_t lo_item_idx = bisect(cumsum, count, idx_lo);
        size_t hi_item_idx = bisect(cumsum, count, idx_hi);

        double x_lo = items[lo_item_idx];
        double x_hi = items[hi_item_idx];
        out[q] = x_lo + frac * (x_hi - x_lo);
    }

    free(cumsum);
    free(items);
    free(weights);
    return 0;
}

static void smooth_log_slope_box(const double *in, size_t n, size_t window_points, double *out)
{
    if (n == 0) {
        return;
    }

    if (window_points <= 1 || n == 1) {
        memcpy(out, in, n * sizeof(double));
        return;
    }

    size_t half = window_points / 2;
    for (size_t i = 0; i < n; i++) {
        size_t start = i > half ? i - half : 0;
        size_t end = i + half;
        if (end >= n) {
            end = n - 1;
        }

        double sum = 0.0;
        for (size_t j = start; j <= end; j++) {
            sum += in[j];
        }

        out[i] = sum / (double)(end - start + 1);
    }
}

int fc_kll_compute_density_from_quantiles(
    const fc_kll_sketch *sketch,
    size_t grid_points,
    fc_kll_density *out,
    bool smooth,
    size_t smoothing_points)
{
    if (sketch == NULL || out == NULL || grid_points < 2
        || (smooth && (smoothing_points == 0 || (smoothing_points % 2) == 0))) {
        return -1;
    }

    out->support = NULL;
    out->values = NULL;
    out->n = 0;

    if (sketch->n == 0) {
        return 0;
    }

    double *ps = fc_alloc(grid_points, sizeof(double), "kll_qgrid_ps");
    double *xs = fc_alloc(grid_points, sizeof(double), "kll_qgrid_xs");
    for (size_t i = 0; i < grid_points; i++) {
        ps[i] = (double)i / (double)(grid_points - 1);
    }

    if (fc_kll_compute_quantiles(sketch, ps, grid_points, xs) != 0) {
        free(ps);
        free(xs);
        return -1;
    }
    free(ps);

    // remove duplicates in the quantile grid to avoid zero-width intervals
    double *ux = fc_alloc(grid_points, sizeof(double), "kll_qgrid_ux");
    double *up = fc_alloc(grid_points, sizeof(double), "kll_qgrid_up");
    size_t unique = 0;
    for (size_t i = 0; i < grid_points; i++) {
        double p = (double)i / (double)(grid_points - 1);
        if (unique == 0) {
            ux[unique] = xs[i];
            up[unique] = p;
            unique = 1;
            continue;
        }

        if (xs[i] > ux[unique - 1]) {
            ux[unique] = xs[i];
            up[unique] = p;
            unique++;
        } else {
            up[unique - 1] = p;
        }
    }
    free(xs);

    if (unique == 0) {
        free(ux);
        free(up);
        return 0;
    }

    if (unique == 1) {
        out->support = fc_alloc(1, sizeof(double), "kll_qgrid_density_support_single");
        out->values = fc_alloc(1, sizeof(double), "kll_qgrid_density_values_single");
        out->support[0] = ux[0];
        out->values[0] = 1.0;
        out->n = 1;
        free(ux);
        free(up);
        return 0;
    }

    // compute density values for each interval between unique quantiles
    // we do smoothing in log-space of the derivative to avoid negative densities
    size_t max_intervals = unique - 1;
    out->support = fc_alloc(max_intervals, sizeof(double), "kll_qgrid_density_support");
    out->values = fc_alloc(max_intervals, sizeof(double), "kll_qgrid_density_values");
    double *widths = fc_alloc(max_intervals, sizeof(double), "kll_qgrid_density_widths");
    double *log_slope = smooth ? fc_alloc(max_intervals, sizeof(double), "kll_qgrid_density_log_slope") : NULL;
    double *log_slope_smoothed = smooth ? fc_alloc(max_intervals, sizeof(double), "kll_qgrid_density_log_slope_smooth") : NULL;

    size_t n_out = 0;
    for (size_t i = 0; i < max_intervals; i++) {
        double dp = up[i + 1] - up[i];
        double dx = ux[i + 1] - ux[i];
        if (!(dp > 0.0) || !(dx > 0.0)) {
            continue;
        }

        out->support[n_out] = 0.5 * (ux[i] + ux[i + 1]);
        widths[n_out] = dx;
        if (smooth) {
            log_slope[n_out] = log(dx / dp);
        } else {
            out->values[n_out] = dp / dx;
        }
        n_out++;
    }

    if (smooth && n_out > 0) {
        // smoothing is just moving average in log-space of the slope (dx/dp)
        // Tried fancier filters, but they weren't visually better for final result.
        // In practice the user may want to smooth end result with a kernel of their choice.
        smooth_log_slope_box(log_slope, n_out, smoothing_points, log_slope_smoothed);
        for (size_t i = 0; i < n_out; i++) {
            out->values[i] = exp(-log_slope_smoothed[i]);
        }
    }

    // normalize density values
    if (n_out > 0) {
        double mass = 0.0;
        for (size_t i = 0; i < n_out; i++) {
            mass += out->values[i] * widths[i];
        }
        if (mass > 0.0) {
            double inv_mass = 1.0 / mass;
            for (size_t i = 0; i < n_out; i++) {
                out->values[i] *= inv_mass;
            }
        }
    }

    out->n = n_out;

    free(log_slope_smoothed);
    free(log_slope);
    free(widths);
    free(ux);
    free(up);
    return 0;
}

void fc_kll_density_destroy(fc_kll_density *density)
{
    if (density == NULL) {
        return;
    }
    free(density->support);
    free(density->values);
    density->support = NULL;
    density->values = NULL;
    density->n = 0;
}
