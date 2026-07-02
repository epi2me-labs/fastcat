#ifndef FASTCAT_KLL_H
#define FASTCAT_KLL_H

#include <stdbool.h>
#include <stddef.h>

// KLL (Karnin-Lang-Liberty) sketch for memory-efficient quantile computation
//
// Memory: O(k log n) where k is accuracy parameter, n is number of items
// Error bound: ~1/k, so k=200 gives 1% rank error
//
// Typical values:
//   k=200: ~2KB memory, 1.0% rank error
//   k=400: ~4KB memory, 0.5% rank error

typedef struct fc_kll_sketch fc_kll_sketch;

// Quantile triple (Q25, Q50, Q75) for distribution summary
typedef struct {
    double q25;
    double q50;
    double q75;
} fc_kll_quantile_triple;

// Density estimate with explicit support/value arrays.
// support[i] is the i-th support location.
// values[i] is the estimated density at support[i].
// n is the number of support/value entries.
typedef struct {
    double *support;
    double *values;
    size_t n;
} fc_kll_density;

// Create a new KLL sketch with accuracy parameter k
// Larger k gives tighter error bounds but uses more memory
// Typical range: k=100 to k=400
fc_kll_sketch *fc_kll_create(size_t k);

// Add a value to the sketch
void fc_kll_add(fc_kll_sketch *sketch, double value);

// Compute one or more quantiles from sketch
// sketch: KLL sketch to query
// ps: array of quantile values (each in [0, 1])
// n: number of quantiles to compute
// out: output array (must be pre-allocated with n elements)
// Returns: 0 on success, -1 on error
int fc_kll_compute_quantiles(const fc_kll_sketch *sketch, const double *ps, size_t n, double *out);

// Convenience wrapper: compute a single quantile (p in [0, 1])
// Returns: interpolated quantile value, or NAN if sketch is empty or p is invalid
double fc_kll_compute_quantile(const fc_kll_sketch *sketch, double p);

// Convenience wrapper: compute Q25/Q50/Q75 in one pass
// Returns: quantile triple, or NAN values if sketch is empty
fc_kll_quantile_triple fc_kll_compute_quantile_triple(const fc_kll_sketch *sketch);

// Compute density from a uniform quantile-probability grid.
// grid_points controls resolution in probability space (must be >= 2).
// If smooth is true, a centered box smoothing window is applied in log-slope
// space with smoothing_points points (must be odd and >= 1) before converting
// back to density. Output is normalized to unit mass.
// Returns 0 on success, -1 on invalid input.
int fc_kll_compute_density_from_quantiles(
    const fc_kll_sketch *sketch,
    size_t grid_points,
    fc_kll_density *out,
    bool smooth,
    size_t smoothing_points);

// Free arrays owned by a density estimate and reset fields.
void fc_kll_density_destroy(fc_kll_density *density);

// Get total number of items added to sketch
size_t fc_kll_count(const fc_kll_sketch *sketch);

// Return bytes allocated by the sketch data structure itself.
// Includes the sketch struct, level array capacity, and per-level item buffers.
// Excludes allocator bookkeeping overhead.
size_t fc_kll_memory_bytes(const fc_kll_sketch *sketch);

// Free sketch memory
void fc_kll_destroy(fc_kll_sketch *sketch);

#endif
