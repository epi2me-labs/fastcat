#include "core/fc_stats.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"

// Bin specification: array of {upper_bound, step_size} pairs
typedef struct {
    size_t upper;
    size_t step;
} bin_segment;

// Create segmented length histogram from bin specification
static read_stats *create_segmented_length_stats(const bin_segment *segments, size_t n_segments, const char *label)
{
    read_stats *stats = fc_alloc(1, sizeof(read_stats), label);

    bin_groups *bins = fc_alloc(1, sizeof(bin_groups), "bin_groups");
    stats->buckets = bins;
    bins->n = n_segments;
    bins->groups = fc_alloc(n_segments, 3 * sizeof(size_t), "groups");

    // Copy segment specification into groups array
    for (size_t i = 0; i < n_segments; i++) {
        bins->groups[3 * i] = segments[i].upper;
        bins->groups[3 * i + 1] = segments[i].step;
        bins->groups[3 * i + 2] = 0; // computed below
    }

    // Count the total number of bins across all groups
    stats->n = 0;
    size_t lower = 0;
    for (size_t i = 0; i < bins->n; i++) {
        size_t upper = bins->groups[3 * i];
        size_t step = bins->groups[3 * i + 1];
        if (step == 0) {
            FC_LOG_FATAL("Invalid histogram configuration: step size cannot be zero");
        }
        size_t nbins = (upper - lower) / step;
        bins->groups[3 * i + 2] = nbins;
        stats->n += nbins;
        lower = upper;
    }
    stats->n++; // overflow bin

    // Fill in all the edges
    stats->width = 0;
    // Keep n+1 edges so the overflow bucket [edges[n-1], edges[n]) is representable.
    stats->edges = fc_alloc(stats->n + 1, sizeof(size_t), "edges");
    stats->counts = fc_alloc(stats->n, sizeof(size_t), "counts");
    size_t i = 0;
    lower = 0;
    size_t upper = 0;
    for (size_t b = 0; b < bins->n; b++) {
        upper = bins->groups[3 * b];
        size_t step = bins->groups[3 * b + 1];
        for (size_t j = lower; j < upper; j += step, i++) {
            stats->edges[i] = j;
        }
        lower = upper;
    }
    stats->edges[i] = upper;
    // Final edge for the overflow bucket upper bound.
    stats->edges[stats->n] = SIZE_MAX;
    return stats;
}

read_stats *create_length_stats(void)
{
    // Single massive bucket up to 10M (76MB per histogram)
    bin_segment spec[] = {
        { 10000000, 1 }
    };
    return create_segmented_length_stats(spec, 1, "length_stats");
}

read_stats *create_qc_length_bins(void)
{
    // QC-optimized segmented binning.
    // The point here is that we may create many of these, one per time window,
    // or qscore bin, so we can't have them all taking 76MB. Ultimately we are
    // just going to be plotting a heatmap/boxplot so we don't need complete
    // resolution.
    // For quantiles, we could replace this with KLL sketch which has
    // provable error bounds.
    // Memory per histogram: 6,351 * 8 bytes = ~51 KB
    bin_segment spec[] = {
        { 5000, 1 }, // precise for short reads
        { 10000, 10 }, // good for medium reads
        { 50000, 50 }, // adequate for fat tail
        { 100000, 1000 } // whateva bruv
    };
    return create_segmented_length_stats(spec, 4, "qc_length_stats");
}

void destroy_length_stats(read_stats *stats)
{
    if (stats != NULL) {
        free(stats->buckets->groups);
        free(stats->buckets);
        free(stats->edges);
        free(stats->counts);
        free(stats);
    }
}

void add_length_count(read_stats *stats, size_t x)
{
    size_t lower = 0;
    size_t cum_bin = 0;
    bool done = false;
    for (size_t i = 0; i < stats->buckets->n; i++) {
        size_t upper = stats->buckets->groups[3 * i];
        if (x < upper) {
            size_t step = stats->buckets->groups[3 * i + 1];
            stats->counts[cum_bin + (x - lower) / step]++;
            done = true;
            break;
        }
        lower = upper;
        cum_bin += stats->buckets->groups[3 * i + 2];
    }
    if (!done) {
        stats->counts[cum_bin]++;
    }
}

read_stats *create_qual_stats(float width)
{
    if (width <= 0.0f || !isfinite(width)) {
        FC_LOG_FATAL("Invalid histogram width: %f (must be positive and finite)", width);
    }
    read_stats *stats = fc_alloc(1, sizeof(read_stats), "quality stats");
    stats->width = width;
    // this fixes the range to [0, 100], good for both QUAL and %age acc
    stats->n = (size_t)(100.0 / stats->width) + 1;

    // Validate resulting size is reasonable (prevent massive allocations)
    if (stats->n > 10000000) {
        FC_LOG_FATAL("Histogram width %f results in excessive bins (%zu)", width, stats->n);
    }

    stats->counts = fc_alloc(stats->n, sizeof(size_t), "counts");
    return stats;
}

void destroy_qual_stats(read_stats *stats)
{
    if (stats != NULL) {
        free(stats->counts);
        free(stats);
    }
}

void add_qual_count(read_stats *stats, float q)
{
    // Non-finite (NaN, Inf), negative, and zero quality scores bin to index 0
    size_t idx = 0;
    if (isfinite(q) && q > 0.0f) {
        if (q >= 100.0f) {
            idx = stats->n - 1;
        } else {
            idx = (size_t)(q / stats->width);
            if (idx >= stats->n) {
                idx = stats->n - 1;
            }
        }
    }
    stats->counts[idx]++;
}

// Counts leading decimal places for histogram bin formatting.
// E.g., 0.001 -> 3, 0.02 -> 2. Uses string formatting to avoid
// floating-point arithmetic edge cases.
size_t _leading_decimals(float number)
{
    char str[64] = { 0 };
    FC_SNPRINTF_CHECK(str, sizeof(str), "%f", number);
    char *point = strchr(str, '.');
    return 1 + strspn(point + 1, "0");
}

void print_stats(read_stats *stats, bool zeroes, bool tsv, FILE *fp)
{
    if (fp == NULL) {
        fp = stderr;
    }
    if (stats->width == 0) {
        for (size_t i = 0; i < stats->n; i++) {
            if (stats->counts[i] == 0 && !zeroes) {
                continue;
            }
            if (tsv) {
                fprintf(fp, "%zu\t%zu\t%zu\n", stats->edges[i], stats->edges[i + 1], stats->counts[i]);
            } else {
                fprintf(fp, "[%zu, %zu)\t%zu\n", stats->edges[i], stats->edges[i + 1], stats->counts[i]);
            }
        }
    } else {
        if (tsv) {
            size_t decimals = _leading_decimals(stats->width);
            char fmt[64] = { 0 }; // my brain hurts and 64 seems big enough
            FC_SNPRINTF_CHECK(fmt, sizeof(fmt), "%%.%zuf\t%%.%zuf\t%%zu\n", decimals, decimals);
            for (size_t i = 0; i < stats->n; i++) {
                if (stats->counts[i] == 0 && !zeroes) {
                    continue;
                }
                fprintf(fp, fmt, (float)i * stats->width, (float)(i + 1) * stats->width, stats->counts[i]);
            }
        } else {
            size_t decimals = _leading_decimals(stats->width);
            char fmt[64] = { 0 };
            FC_SNPRINTF_CHECK(fmt, sizeof(fmt), "[%%.%zuf, %%.%zuf)\t%%zu\n", decimals, decimals);
            for (size_t i = 0; i < stats->n; i++) {
                if (stats->counts[i] == 0 && !zeroes) {
                    continue;
                }
                fprintf(fp, fmt, (float)i * stats->width, (float)(i + 1) * stats->width, stats->counts[i]);
            }
        }
    }
}
