#ifndef FASTCAT_QUANTILE_H
#define FASTCAT_QUANTILE_H

#include <stddef.h>

#include "core/fc_stats.h"

// Quantile triple (Q25, Q50, Q75) for distribution summary
typedef struct {
    double q25;
    double q50;
    double q75;
} fc_quantile_triple;

// Compute one or more quantiles from histogram using linear interpolation over CDF
// stats: histogram to compute quantiles from
// p: array of quantile values (each 0.0 to 1.0, e.g., 0.5 for median)
// n: number of quantiles to compute
// out: output array (must be pre-allocated with n elements)
// Returns: 0 on success, -1 on error
int fc_compute_quantiles(const read_stats *stats, const float *p, size_t n, double *out);

// Convenience wrapper: compute a single quantile
// Returns: interpolated quantile value, or NAN if histogram is empty or invalid
double fc_compute_quantile(const read_stats *stats, float p);

// Convenience wrapper: compute Q25/Q50/Q75 in one pass
// Returns: quantile triple, or NAN values if histogram is empty
fc_quantile_triple fc_compute_quantile_triple(const read_stats *stats);

#endif
