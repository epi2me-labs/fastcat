#include "core/fc_quantile.h"

#include <math.h>
#include <stdlib.h>

#include "core/fc_common.h"

// get bin edges for interpolation
static inline void get_bin_edges(const read_stats *stats, size_t bin_idx, double *lower, double *upper)
{
    if (stats->width > 0) {
        // fixed-width bins
        *lower = (double)bin_idx * stats->width;
        *upper = (double)(bin_idx + 1) * stats->width;
    } else {
        // segmented bins
        *lower = (double)stats->edges[bin_idx];
        *upper = (double)stats->edges[bin_idx + 1];
    }
}

int fc_compute_quantiles(const read_stats *stats, const float *p, size_t n, double *out)
{
    if (!stats || !p || !out || n == 0) {
        return -1;
    }

    // build cumulative sum array
    size_t *cumsum = fc_alloc(stats->n, sizeof(size_t), "cumsum");
    cumsum[0] = stats->counts[0];
    for (size_t i = 1; i < stats->n; i++) {
        cumsum[i] = cumsum[i - 1] + stats->counts[i];
    }
    size_t total_count = cumsum[stats->n - 1];

    // handle empty histogram
    if (total_count == 0) {
        for (size_t i = 0; i < n; i++) {
            out[i] = NAN;
        }
        free(cumsum);
        return 0;
    }

    // compute each quantile
    for (size_t q = 0; q < n; q++) {
        // validate quantile
        if (p[q] < 0.0f || p[q] > 1.0f) {
            out[q] = NAN;
            continue;
        }

        double target = p[q] * (double)total_count;

        // p == 0 should resolve to the lower edge of the first non-empty bin.
        // The binary-search predicate below assumes a strictly positive target
        // and otherwise can stop on an empty leading bucket.
        if (target <= 0.0) {
            size_t bin_idx = 0;
            while (stats->counts[bin_idx] == 0) {
                bin_idx++;
            }
            double lower, upper;
            get_bin_edges(stats, bin_idx, &lower, &upper);
            (void)upper;
            out[q] = lower;
            continue;
        }

        // binary search for bin containing target
        size_t left = 0, right = stats->n - 1;
        size_t bin_idx = 0;

        while (left <= right) {
            size_t mid = left + (right - left) / 2;
            size_t cumsum_before = (mid == 0) ? 0 : cumsum[mid - 1];

            if ((double)cumsum_before < target && (double)cumsum[mid] >= target) {
                bin_idx = mid;
                break;
            } else if ((double)cumsum[mid] < target) {
                left = mid + 1;
            } else {
                if (mid == 0) {
                    break;
                }
                right = mid - 1;
            }
        }

        // we've found the bin, we could just take the lower edge, but instead
        // we can assume a uniform distribution of values within the bin and
        // interpolate to get a better estimate.
        double lower, upper;
        get_bin_edges(stats, bin_idx, &lower, &upper);

        // The segmented overflow bucket and the terminal fixed-width bucket are
        // both open-ended in practice, so return the lower bound rather than
        // inventing values beyond the histogram's representable range.
        if (stats->width == 0 && stats->edges[bin_idx + 1] == SIZE_MAX) {
            out[q] = lower;
            continue;
        }
        if (stats->width > 0 && bin_idx == stats->n - 1) {
            out[q] = lower;
            continue;
        }

        // how far through the bin is our target quantile?
        size_t bin_count = stats->counts[bin_idx];
        size_t cumsum_before = (bin_idx == 0) ? 0 : cumsum[bin_idx - 1];
        double pos_in_bin = (target - (double)cumsum_before) / (double)bin_count;
        out[q] = lower + pos_in_bin * (upper - lower);
    }

    free(cumsum);
    return 0;
}

// convenience wrapper for single quantile
double fc_compute_quantile(const read_stats *stats, float p)
{
    double result;
    if (fc_compute_quantiles(stats, &p, 1, &result) != 0) {
        return NAN;
    }
    return result;
}

// convenience wrapper for Q25/Q50/Q75
fc_quantile_triple fc_compute_quantile_triple(const read_stats *stats)
{
    float p[3] = { 0.25f, 0.50f, 0.75f };
    double out[3];
    fc_quantile_triple result;

    if (fc_compute_quantiles(stats, p, 3, out) != 0) {
        result.q25 = result.q50 = result.q75 = NAN;
    } else {
        result.q25 = out[0];
        result.q50 = out[1];
        result.q75 = out[2];
    }

    return result;
}
