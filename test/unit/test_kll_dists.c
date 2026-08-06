/*
 * Stress tests for KLL quantile estimation across various distributions.
 *
 * These tests generate large samples from known distributions, add them to a KLL sketch,
 * and then compute quantiles. Results are validated using both rank error and value error.
 *
 * 1. Rank Error, distribution independent (tests KLL guarantee):
 *    - Rank error = |r - p|
 *    - Assert: rank_error < 1.5% (KLL with k=200 guarantees error ~ 1%)
 *
 * 2. Value Error, distribution dependent (tests practical accuracy):
 *    - Value error = |x - v_exact|
 *    - Assert: value_error < adaptive tolerance (2% base, up to 3.5x at extremes)
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef _WIN32
#include <sys/resource.h>
#endif

#include "unity.h"
#include "main.h"
#include "core/fc_kll.h"

// Global profiling state
static double total_kll_time_sec = 0.0;
static size_t peak_kll_memory_bytes = 0;
static size_t baseline_memory_bytes = 0;
static size_t total_sketches_created = 0;
static size_t total_values_added = 0;
static size_t total_quantiles_computed = 0;

static size_t get_current_rss_bytes(void)
{
#ifdef _WIN32
    /* The stress test's assertions do not depend on process RSS. */
    return 0;
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    // ru_maxrss is in bytes on macOS, kilobytes on Linux
#ifdef __APPLE__
    return (size_t)usage.ru_maxrss;
#else
    return (size_t)usage.ru_maxrss * 1024;
#endif
#endif
}

typedef struct {
    const char *name;
    double (*quantile_fn)(double p, double param1, double param2); // inverse CDF
    double (*cdf_fn)(double x, double param1, double param2); // CDF for rank error testing
    double param1;
    double param2; // unused in some distributions
} distribution_fixture;

// Uniform(min, max)
static double uniform_quantile(double p, double min, double max)
{
    return min + p * (max - min);
}

static double uniform_cdf(double x, double min, double max)
{
    if (x <= min) return 0.0;
    if (x >= max) return 1.0;
    return (x - min) / (max - min);
}

// Exponential(lambda)
static double exponential_quantile(double p, double lambda, double unused)
{
    (void)unused;
    return -log(1.0 - p) / lambda;
}

static double exponential_cdf(double x, double lambda, double unused)
{
    (void)unused;
    if (x < 0.0) return 0.0;
    return 1.0 - exp(-lambda * x);
}

// Logistic(mu, s)
static double logistic_quantile(double p, double mu, double s)
{
    return mu + s * log(p / (1.0 - p));
}

static double logistic_cdf(double x, double mu, double s)
{
    return 1.0 / (1.0 + exp(-(x - mu) / s));
}

// Cauchy(x0, gamma)
static double cauchy_quantile(double p, double x0, double gamma)
{
    return x0 + gamma * tan(M_PI * (p - 0.5));
}

static double cauchy_cdf(double x, double x0, double gamma)
{
    return 0.5 + atan((x - x0) / gamma) / M_PI;
}

// Weibull(lambda, k)
static double weibull_quantile(double p, double lambda, double k)
{
    return lambda * pow(-log(1.0 - p), 1.0 / k);
}

static double weibull_cdf(double x, double lambda, double k)
{
    if (x < 0.0) return 0.0;
    return 1.0 - exp(-pow(x / lambda, k));
}

// Pareto(xm, alpha)
static double pareto_quantile(double p, double xm, double alpha)
{
    return xm / pow(1.0 - p, 1.0 / alpha);
}

static double pareto_cdf(double x, double xm, double alpha)
{
    if (x < xm) return 0.0;
    return 1.0 - pow(xm / x, alpha);
}

// Laplace(mu, b)
static double laplace_quantile(double p, double mu, double b)
{
    if (p < 0.5) {
        return mu + b * log(2.0 * p);
    } else {
        return mu - b * log(2.0 * (1.0 - p));
    }
}

static double laplace_cdf(double x, double mu, double b)
{
    if (x < mu) {
        return 0.5 * exp((x - mu) / b);
    } else {
        return 1.0 - 0.5 * exp(-(x - mu) / b);
    }
}

// simple deterministic RNG for reproducibility
static uint32_t rng_state = 42;

static uint32_t xorshift32(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

// Fisher-Yates shuffle
static void shuffle_doubles(double *arr, size_t n)
{
    for (size_t i = n - 1; i > 0; i--) {
        size_t j = xorshift32() % (i + 1);
        double tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

static void run_distribution_quantiles(distribution_fixture dist, size_t n_samples)
{
    rng_state = 42;

    // generate exact samples from distribution using inverse CDF.
    // This uniformlly samples in probability space, then transforms
    // to value space using the quantile function (inverse CDF).
    // The same function is used later to check. What goes in should
    // come out. Roll Safe.
    double *samples = malloc(n_samples * sizeof(double));
    for (size_t i = 1; i < n_samples; i++) {
        double p = (double)i / (double)n_samples;
        samples[i - 1] = dist.quantile_fn(p, dist.param1, dist.param2);
    }

    // shuffle to simulate realistic unsorted data
    shuffle_doubles(samples, n_samples - 1);

    // add to KLL sketch - profile time and memory
    clock_t start = clock();
    size_t mem_before = get_current_rss_bytes();

    fc_kll_sketch *sketch = fc_kll_create(200);
    total_sketches_created++;

    for (size_t i = 0; i < n_samples - 1; i++) {
        fc_kll_add(sketch, samples[i]);
        total_values_added++;
    }
    free(samples);

    // Measure actual memory growth
    size_t mem_after = get_current_rss_bytes();
    if (mem_after > peak_kll_memory_bytes) {
        peak_kll_memory_bytes = mem_after;
    }
    (void)mem_before; // suppress unused warning - we only care about peak

    // test quantiles at various points
    // Skip extreme tails (p<0.05, p>0.95) where sample density is low
    double test_ps[] = {0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95};
    size_t n_test = sizeof(test_ps) / sizeof(test_ps[0]);
    double kll_out[9];

    fc_kll_compute_quantiles(sketch, test_ps, n_test, kll_out);
    total_quantiles_computed += n_test;

    clock_t end = clock();
    total_kll_time_sec += (double)(end - start) / CLOCKS_PER_SEC;

    // Validate using BOTH rank error and value error
    for (size_t i = 0; i < n_test; i++) {
        double target_p = test_ps[i];
        double exact_value = dist.quantile_fn(target_p, dist.param1, dist.param2);
        double kll_estimate = kll_out[i];

        // Rank error (distribution-independent) - 1%
        double actual_rank = dist.cdf_fn(kll_estimate, dist.param1, dist.param2);
        double rank_error = fabs(actual_rank - target_p);
        double rank_tolerance = 0.01;

        if (rank_error > rank_tolerance) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                "%s n=%zu: RANK ERROR - p=%.2f, kll=%.3f, actual_rank=%.4f, "
                "rank_error=%.4f (%.2f%%), tolerance=%.2f%%",
                dist.name, n_samples, target_p, kll_estimate, actual_rank,
                rank_error, rank_error * 100.0, rank_tolerance * 100.0);
            TEST_FAIL_MESSAGE(msg);
        }

        // Value error (distribution-dependent, tests practical accuracy)
        double value_error = fabs(kll_estimate - exact_value);
        double extremeness = fabs(target_p - 0.5) * 2.0; // 0 at median, 1 at extremes
        double base_tolerance = fmax(fabs(exact_value) * 0.02, 1.2); // 2% or 1.2
        double value_tolerance = base_tolerance * (1.0 + extremeness * 2.5); // up to 3.5x

        if (value_error > value_tolerance) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                "%s n=%zu: VALUE ERROR - p=%.2f, exact=%.3f, kll=%.3f, "
                "value_error=%.3f, tolerance=%.3f",
                dist.name, n_samples, target_p, exact_value, kll_estimate,
                value_error, value_tolerance);
            TEST_FAIL_MESSAGE(msg);
        }
    }

    fc_kll_destroy(sketch);
}

static void run_dist_quantiles_fixture_list(
    const distribution_fixture *dists, size_t n_dists)
{
    for (size_t i = 0; i < n_dists; i++) {
        for (size_t n = 1000; n <= 1000000; n *= 10) {
            run_distribution_quantiles(dists[i], n);
        }
    }
}

// In each case below the parameter choices are simply those used on the
// Wikipedia page for that distribution in the example plots. They are not
// chosen for any other reason.
#define DIST_FIXTURE_1(name, quantile_fn, cdf_fn, p1_label, p1) \
    { name "(" p1_label "=" #p1 ")", quantile_fn, cdf_fn, (p1), 0.0 }

#define DIST_FIXTURE_2(name, quantile_fn, cdf_fn, p1_label, p1, p2_label, p2) \
    { name "(" p1_label "=" #p1 "," p2_label "=" #p2 ")", quantile_fn, cdf_fn, (p1), (p2) }

#define Uniform(a, b) \
    { "Uniform(" #a "," #b ")", uniform_quantile, uniform_cdf, (a), (b) }
static void uniform_stress(void)
{
    distribution_fixture dists[] = {
        // ok, well this one is silly
        Uniform(0, 1),
        Uniform(5, 10),
        Uniform(0, 100),
        Uniform(-100, 1000),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Exponential(lambda) \
    DIST_FIXTURE_1("Exponential", exponential_quantile, exponential_cdf, "lambda", lambda)
static void exponential_stress(void)
{
    distribution_fixture dists[] = {
        Exponential(0.1),
        Exponential(0.5),
        Exponential(1.0),
        Exponential(2.0),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Logistic(mu, s) \
    DIST_FIXTURE_2("Logistic", logistic_quantile, logistic_cdf, "mu", mu, "s", s)
static void logistic_stress(void)
{
    distribution_fixture dists[] = {
        Logistic(5, 2),
        Logistic(9, 3),
        Logistic(9, 4),
        Logistic(6, 2),
        Logistic(2, 1),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Weibull(lambda, k) \
    DIST_FIXTURE_2("Weibull", weibull_quantile, weibull_cdf, "lambda", lambda, "k", k)
static void weibull_stress(void)
{
    distribution_fixture dists[] = {
        Weibull(1, 0.5),
        Weibull(1, 1),
        Weibull(1, 1.5),
        Weibull(1, 5),
        // Wiki just uses lambda=1, so lets also do something different
        Weibull(10, 0.5),
        Weibull(10, 1),
        Weibull(10, 1.5),
        Weibull(10, 5),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Cauchy(x0, gamma) \
    DIST_FIXTURE_2("Cauchy", cauchy_quantile, cauchy_cdf, "x0", x0, "gamma", gamma)
static void cauchy_stress(void)
{
    distribution_fixture dists[] = {
        Cauchy(0, 0.5),
        Cauchy(0, 1),
        Cauchy(0, 2),
        Cauchy(-2, 1),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Pareto(xm, alpha) \
    DIST_FIXTURE_2("Pareto", pareto_quantile, pareto_cdf, "xm", xm, "alpha", alpha)
static void pareto_stress(void)
{
    distribution_fixture dists[] = {
        Pareto(1, 1),
        Pareto(1, 2),
        Pareto(1, 3),
        Pareto(1, 9999999),  // Good ol' Dirac
        // Wiki just uses xm=1, so lets also do something different
        Pareto(2, 1),
        Pareto(2, 2),
        Pareto(2, 3),
        Pareto(2, 9999999),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

#define Laplace(mu, b) \
    DIST_FIXTURE_2("Laplace", laplace_quantile, laplace_cdf, "mu", mu, "b", b)
static void laplace_stress(void)
{
    distribution_fixture dists[] = {
        Laplace(0, 1),
        Laplace(0, 2),
        Laplace(0, 4),
        Laplace(-5, 4),
    };
    run_dist_quantiles_fixture_list(dists, sizeof(dists) / sizeof(dists[0]));
}

FC_TEST(quantile_stress_suite, test_quantile_stress_battery)
{
    // reset profiling state at start
    total_kll_time_sec = 0.0;
    peak_kll_memory_bytes = 0;
    baseline_memory_bytes = get_current_rss_bytes();
    total_sketches_created = 0;
    total_values_added = 0;
    total_quantiles_computed = 0;

    uniform_stress();
    exponential_stress();
    logistic_stress();
    weibull_stress();
    cauchy_stress();
    pareto_stress();
    laplace_stress();

    // print profiling results
    size_t mem_growth = peak_kll_memory_bytes - baseline_memory_bytes;
    fprintf(stderr, "========== KLL Sketch Profiling ==========\n");
    fprintf(stderr, "Total sketches created:     %zu\n", total_sketches_created);
    fprintf(stderr, "Total values added:         %.2f M\n", total_values_added / 1e6);
    fprintf(stderr, "Total quantiles computed:   %.2f K\n", total_quantiles_computed / 1e3);
    fprintf(stderr, "Memory growth (peak-base):  %.1f KB\n", mem_growth / 1024.0);
    fprintf(stderr, "Total CPU time:             %.3f sec\n", total_kll_time_sec);
    fprintf(stderr, "Throughput (add):           %.2f M values/sec\n",
        total_values_added / 1e6 / total_kll_time_sec);
    fprintf(stderr, "Throughput (quantile):      %.2f K queries/sec\n",
        total_quantiles_computed / 1e3 / total_kll_time_sec);
    fprintf(stderr, "==========================================\n");
}
