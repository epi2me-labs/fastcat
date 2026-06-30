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

// compare function for qsort
static int compare_doubles(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) {
        return -1;
    }
    if (da > db) {
        return 1;
    }
    return 0;
}

// temporary struct for sorting items with their weights
typedef struct {
    double item;
    size_t weight;
} item_weight_pair;

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
    qsort(comp->items, comp->count, sizeof(double), compare_doubles);
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

size_t fc_kll_count(const fc_kll_sketch *sketch)
{
    return sketch == NULL ? 0 : sketch->n;
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

fc_kll_quantile_triple fc_kll_compute_quantile_triple(const fc_kll_sketch *sketch)
{
    double p[3] = { 0.25, 0.50, 0.75 };
    double out[3] = { NAN, NAN, NAN };
    fc_kll_quantile_triple result;

    fc_kll_compute_quantiles(sketch, p, 3, out);

    result.q25 = out[0];
    result.q50 = out[1];
    result.q75 = out[2];

    return result;
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

    // fetch each quantile through linear search,
    // we could do bisection, or assume inputs are sorted but this
    // is simple and fast enough for small numbers
    for (size_t q = 0; q < n; q++) {
        if (ps[q] < 0.0 || ps[q] > 1.0) {
            out[q] = NAN;
            continue;
        }
        double target_rank = ps[q] * (double)total_weight;
        size_t i;
        for (i = 0; i < count; i++) {
            if ((double)cumsum[i] >= target_rank) {
                out[q] = items[i];
                break;
            }
        }
        if (i == count) {
            out[q] = items[count - 1];
        }
    }

    free(cumsum);
    free(items);
    free(weights);
    return 0;
}
