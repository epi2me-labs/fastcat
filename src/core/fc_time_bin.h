#ifndef FASTCAT_FC_TIME_BIN_H
#define FASTCAT_FC_TIME_BIN_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_kll.h"
#include "core/fc_read.h"

#define FC_TIME_BIN_DEFAULT_SIZE_S 3600.0
#define FC_TIME_BIN_DEFAULT_K 200

typedef struct fc_time_bin fc_time_bin;

typedef struct {
    double bin_size_s;
    size_t kll_k;
} fc_time_bin_config;

typedef fc_kll_quantile_septuple fc_time_bin_quantiles;

typedef struct {
    double time_bin_start_s;
    double time_bin_end_s;
    size_t read_count;
    size_t base_count;
    size_t cumulative_read_count;
    size_t cumulative_base_count;
    double length_mean;
    fc_time_bin_quantiles length_quantiles;
    double qscore_mean;
    fc_time_bin_quantiles qscore_quantiles;
    double speed_mean;
    fc_time_bin_quantiles speed_quantiles;
    size_t accuracy_read_count;
    double accuracy_mean;
    fc_time_bin_quantiles accuracy_quantiles;
} fc_time_bin_row;

typedef struct {
    fc_time_bin_row *rows;
    size_t n;
} fc_time_bin_table;

typedef struct {
    size_t missing;
    size_t invalid;
    size_t duration_missing;
    size_t duration_invalid;
} fc_time_bin_skip_counts;

fc_time_bin *fc_time_bin_create(size_t barcode_count, double bin_size_s, size_t kll_k);
void fc_time_bin_destroy(fc_time_bin *model);

bool fc_time_bin_parse_timestamp(const char *text, double *out);

void fc_time_bin_add_read(fc_time_bin *model, size_t barcode, const fc_read *read);
void fc_time_bin_add_accuracy(fc_time_bin *model, size_t barcode, const fc_read *read, float accuracy);

bool fc_time_bin_seen(const fc_time_bin *model, size_t barcode);
bool fc_time_bin_accuracy_seen(const fc_time_bin *model, size_t barcode);
fc_time_bin_skip_counts fc_time_bin_skipped_counts(const fc_time_bin *model);

int fc_time_bin_materialize(const fc_time_bin *model, size_t barcode, fc_time_bin_table *out);
int fc_time_bin_materialize_merged(
    const fc_time_bin *model,
    size_t barcode,
    size_t max_rows,
    fc_time_bin_table *out);
void fc_time_bin_table_destroy(fc_time_bin_table *table);

#endif
