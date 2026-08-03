#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_channel_map.h"
#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_iter.h"
#include "core/fc_json.h"
#include "core/fc_kll.h"
#include "core/fc_log.h"
#include "core/fc_metric_qscore_2d.h"
#include "core/fc_sinks.h"
#include "core/kh_counter.h"
#include "sink/common.h"

#define FC_SUMMARY_JSON_SCHEMA_VERSION 1
#define FC_SUMMARY_JSON_SCHEMA_ID "https://github.com/epi2me-labs/fastcat/schema/summary_json.schema.json"

typedef struct {
    kh_counter_t *run_ids;
    kh_counter_t *basecallers;
    size_t read_count;
    size_t base_count;
    size_t min_length;
    size_t max_length;
    double quality_sum;
    size_t quality_count;
    bool seen;
} summary_file_accum;

typedef struct {
    char *path;
    kh_counter_t *run_ids;
    kh_counter_t *basecallers;
    size_t read_count;
    size_t base_count;
    size_t min_length;
    size_t max_length;
    double mean_quality;
} summary_file_entry;

typedef struct {
    bool seen;
    kh_counter_t *run_ids;
    kh_counter_t *basecallers;
    summary_file_entry *files;
    size_t n_files;
    size_t cap_files;
    summary_file_accum current;
} summary_json_barcode;

typedef struct {
    fc_summary_json_sink_opts cfg;
    char *output_root;
    char *sample;
    summary_json_barcode *barcodes;
} summary_json_sink;

static void summary_counter_destroy(kh_counter_t *counter)
{
    kh_counter_destroy(counter);
}

static void summary_file_accum_reset(summary_file_accum *accum)
{
    if (NULL == accum) {
        return;
    }
    summary_counter_destroy(accum->run_ids);
    summary_counter_destroy(accum->basecallers);
    memset(accum, 0, sizeof(*accum));
}

static void summary_file_entry_destroy(summary_file_entry *entry)
{
    if (NULL == entry) {
        return;
    }
    free(entry->path);
    summary_counter_destroy(entry->run_ids);
    summary_counter_destroy(entry->basecallers);
    memset(entry, 0, sizeof(*entry));
}

static void summary_counter_add(kh_counter_t **counter, const char *value)
{
    if (NULL == *counter) {
        *counter = kh_counter_init();
    }
    kh_counter_increment(*counter, (char *)(value ? value : ""));
}

static void summary_file_accum_add_read(summary_file_accum *accum, const fc_read *read)
{
    size_t length = read->bam->core.l_qseq;
    accum->read_count++;
    accum->base_count += length;
    accum->quality_sum += read->meta->mean_q;
    accum->quality_count++;
    if (!accum->seen || length < accum->min_length) {
        accum->min_length = length;
    }
    if (!accum->seen || length > accum->max_length) {
        accum->max_length = length;
    }
    accum->seen = true;
    summary_counter_add(&accum->run_ids, read->meta->runid);
    summary_counter_add(&accum->basecallers, read->meta->basecaller);
}

static void summary_barcode_add_read(summary_json_barcode *state, const fc_read *read)
{
    state->seen = true;
    summary_counter_add(&state->run_ids, read->meta->runid);
    summary_counter_add(&state->basecallers, read->meta->basecaller);
    summary_file_accum_add_read(&state->current, read);
}

static void summary_barcode_append_file(
    summary_json_barcode *state,
    const char *path,
    const fc_file_stats *stats,
    bool use_global_stats)
{
    if (!state->current.seen && !use_global_stats) {
        return;
    }

    if (state->n_files == state->cap_files) {
        size_t new_cap = state->cap_files == 0 ? 4 : state->cap_files * 2;
        state->files = fc_realloc(state->files, new_cap * sizeof(*state->files), "summary_json file entries");
        state->cap_files = new_cap;
    }

    summary_file_entry *entry = &state->files[state->n_files++];
    memset(entry, 0, sizeof(*entry));
    entry->path = fc_strdup(path ? path : "", "summary_json file path");
    entry->run_ids = state->current.run_ids;
    entry->basecallers = state->current.basecallers;
    state->current.run_ids = NULL;
    state->current.basecallers = NULL;

    if (use_global_stats && NULL != stats) {
        entry->read_count = stats->n_seqs;
        entry->base_count = stats->n_bases;
        entry->min_length = stats->n_seqs > 0 ? stats->min_length : 0;
        entry->max_length = stats->n_seqs > 0 ? stats->max_length : 0;
        entry->mean_quality = stats->n_seqs > 0 ? stats->mean_quality : 0.0;
    } else {
        entry->read_count = state->current.read_count;
        entry->base_count = state->current.base_count;
        entry->min_length = state->current.seen ? state->current.min_length : 0;
        entry->max_length = state->current.seen ? state->current.max_length : 0;
        entry->mean_quality = state->current.quality_count > 0
            ? state->current.quality_sum / (double)state->current.quality_count
            : 0.0;
    }

    state->current.read_count = 0;
    state->current.base_count = 0;
    state->current.min_length = 0;
    state->current.max_length = 0;
    state->current.quality_sum = 0.0;
    state->current.quality_count = 0;
    state->current.seen = false;
}

static void summary_barcode_destroy(summary_json_barcode *state)
{
    if (NULL == state) {
        return;
    }
    summary_counter_destroy(state->run_ids);
    summary_counter_destroy(state->basecallers);
    for (size_t i = 0; i < state->n_files; i++) {
        summary_file_entry_destroy(&state->files[i]);
    }
    free(state->files);
    summary_file_accum_reset(&state->current);
    memset(state, 0, sizeof(*state));
}

static const read_stats *summary_read_length_stats(const fc_stats_engine *engine)
{
    return engine ? engine->length_stats : NULL;
}

static const read_stats *summary_quality_stats(const fc_stats_engine *engine)
{
    return engine ? engine->qual_stats : NULL;
}

static const char *summary_channel_device_name(size_t highest_channel)
{
    switch (fc_channel_map_device_for_channel(highest_channel)) {
    case FC_CHANNEL_MAP_FLONGLE:
        return "flongle";
    case FC_CHANNEL_MAP_MINION:
        return "minion";
    case FC_CHANNEL_MAP_PROMETHION:
        return "promethion";
    }
    return "unknown";
}

static void summary_json_set_string_or_null(fc_json_object_t *object, const char *key, const char *value)
{
    if (NULL == value || '\0' == value[0]) {
        fc_json_object_set_null(object, key);
        return;
    }
    fc_json_object_set_string(object, key, value);
}

static void summary_json_set_semicolon_list_or_null(fc_json_object_t *object, const char *key, const char *value)
{
    if (NULL == value || '\0' == value[0]) {
        fc_json_object_set_null(object, key);
        return;
    }

    fc_json_array_t *array = fc_json_object_set_array(object, key);
    const char *item = value;
    while ('\0' != *item) {
        const char *sep = strchr(item, ';');
        size_t len = NULL != sep ? (size_t)(sep - item) : strlen(item);
        if (len > 0) {
            char *copy = fc_strndup(item, len, "summary_json semicolon list");
            fc_json_array_append_string(array, copy);
            free(copy);
        }
        if (NULL == sep) {
            break;
        }
        item = sep + 1;
    }
}

static void summary_json_counter_array(fc_json_object_t *object, const char *key, const kh_counter_t *counter)
{
    if (NULL == counter || 0 == kh_size(counter)) {
        fc_json_object_set_null(object, key);
        return;
    }

    fc_json_array_t *rows = fc_json_object_set_array(object, key);
    for (khiter_t k = 0; k < kh_end(counter); ++k) {
        if (!kh_exist(counter, k)) {
            continue;
        }
        fc_json_value_t *row_value = fc_json_object_value();
        fc_json_object_t *row = fc_json_get_object(row_value);
        fc_json_object_set_string(row, "value", kh_key(counter, k));
        fc_json_object_set_number(row, "count", (double)kh_val(counter, k));
        fc_json_array_append_value(rows, row_value);
    }
}

static void summary_json_files(fc_json_object_t *root, const summary_json_barcode *state)
{
    fc_json_array_t *files = fc_json_object_set_array(root, "files");
    for (size_t i = 0; i < state->n_files; i++) {
        const summary_file_entry *entry = &state->files[i];
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *row = fc_json_get_object(value);
        fc_json_object_set_string(row, "path", entry->path ? entry->path : "");
        fc_json_object_set_number(row, "accepted_read_count", (double)entry->read_count);
        fc_json_object_set_number(row, "accepted_base_count", (double)entry->base_count);
        fc_json_object_set_number(row, "min_length_bp", (double)entry->min_length);
        fc_json_object_set_number(row, "max_length_bp", (double)entry->max_length);
        fc_json_object_set_number(row, "mean_quality", entry->mean_quality);
        summary_json_counter_array(row, "run_id_counts", entry->run_ids);
        summary_json_counter_array(row, "basecaller_counts", entry->basecallers);
        fc_json_array_append_value(files, value);
    }
}

static void summary_json_counts(fc_json_object_t *root, const fc_run_summary_row *row, const summary_json_barcode *state)
{
    fc_json_object_t *counts = fc_json_object_set_object(root, "counts");
    fc_json_object_set_number(counts, "read_count", (double)row->read_count);
    fc_json_object_set_number(counts, "base_count", (double)row->base_count);
    summary_json_counter_array(counts, "run_ids", state->run_ids);
    summary_json_counter_array(counts, "basecallers", state->basecallers);
}

// Root read-length output includes aggregate totals and the yield trace; alignment
// distributions reuse this serializer with zero totals and both optional fields omitted.
static fc_json_value_t *summary_json_length_distribution_value(
    const read_stats *stats,
    size_t total_base_count,
    size_t total_read_count,
    bool include_yield,
    bool include_total_base_count)
{
    if (NULL == stats) {
        return NULL;
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *hist = fc_json_get_object(value);

    // Collect non-zero bins
    size_t *starts_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length starts");
    size_t *ends_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length ends");
    size_t *counts_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length counts");
    size_t *yield_arr = include_yield
        ? fc_alloc(stats->n, sizeof(size_t), "summary_json length yields")
        : NULL;
    size_t n_bins = 0;

    double total = 0.0;
    size_t overflow_count = 0;
    size_t overflow_start = 0;
    if (0.0f == stats->width && NULL != stats->edges) {
        for (size_t i = 0; i < stats->n; i++) {
            size_t count = stats->counts[i];
            total += (double)count;
            if (0 == count) {
                continue;
            }
            size_t start = stats->edges[i];
            size_t end = stats->edges[i + 1];
            if (SIZE_MAX == end) {
                overflow_count = count;
                overflow_start = start;
                continue;
            }
            starts_arr[n_bins] = start;
            ends_arr[n_bins] = end;
            counts_arr[n_bins] = count;
            n_bins++;
        }
    }

    if (include_yield) {
        size_t cumulative_yield = 0;
        for (size_t i = n_bins; i > 0; i--) {
            size_t index = i - 1;
            if (starts_arr[index] != 0 && counts_arr[index] > SIZE_MAX / starts_arr[index]) {
                FC_LOG_FATAL("Length yield exceeds summary.json integer range.");
            }
            size_t bin_yield = starts_arr[index] * counts_arr[index];
            if (bin_yield > SIZE_MAX - cumulative_yield) {
                FC_LOG_FATAL("Cumulative length yield exceeds summary.json integer range.");
            }
            cumulative_yield += bin_yield;
            yield_arr[index] = cumulative_yield;
        }
    }

    // Keep histogram statistics separate from the self-describing table payload.
    fc_json_object_t *data = fc_json_object_set_object(hist, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(data, n_bins, include_yield ? 4 : 3);
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("start", starts_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("end", ends_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("count", counts_arr, n_bins));
    if (include_yield) {
        fc_json_array_append_value(
            columns, fc_json_column_uint_adaptive("yield_at_or_above_bp", yield_arr, n_bins));
    }

    free(starts_arr);
    free(ends_arr);
    free(counts_arr);
    free(yield_arr);

    fc_json_object_set_number(hist, "total_count", total);
    if (include_total_base_count) {
        fc_json_object_set_number(hist, "total_base_count", (double)total_base_count);
        fc_json_object_set_number(
            hist, "mean", total_read_count > 0 ? (double)total_base_count / total_read_count : 0.0);
    }
    fc_json_object_set_number(hist, "median", fc_stats_median_start(stats));
    if (overflow_count > 0) {
        fc_json_object_set_number(hist, "overflow_start_bp", (double)overflow_start);
        fc_json_object_set_number(hist, "overflow_count", (double)overflow_count);
    } else {
        fc_json_object_set_null(hist, "overflow_start_bp");
        fc_json_object_set_null(hist, "overflow_count");
    }
    return value;
}

static fc_json_value_t *summary_json_fixed_distribution_value(const read_stats *stats)
{
    if (NULL == stats) {
        return NULL;
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *hist = fc_json_get_object(value);

    // Collect non-zero bins
    double *starts_arr = fc_alloc(stats->n, sizeof(double), "quality starts");
    double *ends_arr = fc_alloc(stats->n, sizeof(double), "quality ends");
    size_t *counts_arr = fc_alloc(stats->n, sizeof(size_t), "quality counts");
    size_t n_bins = 0;

    double total = 0.0;
    for (size_t i = 0; i < stats->n; i++) {
        size_t count = stats->counts[i];
        total += (double)count;
        if (count == 0) {
            continue;
        }
        starts_arr[n_bins] = i * stats->width;
        ends_arr[n_bins] = (i + 1) * stats->width;
        counts_arr[n_bins] = count;
        n_bins++;
    }

    fc_json_object_set_number(hist, "total_count", total);
    fc_json_object_set_number(hist, "mean", fc_stats_mean_midpoint(stats));
    fc_json_object_set_number(hist, "median", fc_stats_median_start(stats));

    // Keep histogram statistics separate from the self-describing table payload.
    fc_json_object_t *data = fc_json_object_set_object(hist, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(data, n_bins, 3);
    fc_json_array_append_value(columns, fc_json_column_f64("start", starts_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_f64("end", ends_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("count", counts_arr, n_bins));

    free(starts_arr);
    free(ends_arr);
    free(counts_arr);

    return value;
}

static fc_json_value_t *summary_json_heatmap_ndarray(const double *counts, size_t rows, size_t cols)
{
    return fc_json_ndarray_f64_to_uint_adaptive(counts, rows, cols);
}

static void summary_json_length_heatmap(
    fc_json_object_t *heatmaps,
    const fc_length_qscore_2d *model,
    size_t barcode)
{
    if (NULL == model || !fc_length_qscore_2d_seen(model, barcode)) {
        fc_json_object_set_null(heatmaps, "read_length_quality");
        return;
    }

    fc_length_qscore_2d_matrix matrix = { 0 };
    if (fc_length_qscore_2d_materialize(model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize length heatmap for summary.json.");
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *obj = fc_json_get_object(value);

    // Generic 2D array with labeled axes
    fc_json_object_set_string(obj, "x_label", "read_length");
    fc_json_array_t *x_range = fc_json_object_set_array(obj, "x_range");
    fc_json_array_append_number(x_range, matrix.x_lower);
    fc_json_array_append_number(x_range, matrix.x_upper);

    fc_json_object_set_string(obj, "y_label", "qscore");
    fc_json_array_t *y_range = fc_json_object_set_array(obj, "y_range");
    fc_json_array_append_number(y_range, 0.0);
    fc_json_array_append_number(y_range, matrix.qscore_bins * matrix.qscore_width);
    fc_json_object_set_value(
        obj, "data", summary_json_heatmap_ndarray(matrix.counts, matrix.qscore_bins, matrix.x_bins));

    fc_json_object_set_value(heatmaps, "read_length_quality", value);
    fc_length_qscore_2d_matrix_destroy(&matrix);
}

static void summary_json_speed_heatmap(
    fc_json_object_t *heatmaps,
    const fc_speed_qscore_2d *model,
    size_t barcode)
{
    if (NULL == model || !fc_speed_qscore_2d_seen(model, barcode)) {
        fc_json_object_set_null(heatmaps, "translocation_speed_quality");
        return;
    }

    fc_speed_qscore_2d_matrix matrix = { 0 };
    if (fc_speed_qscore_2d_materialize(model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize speed heatmap for summary.json.");
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *obj = fc_json_get_object(value);

    // Generic 2D array with labeled axes
    fc_json_object_set_string(obj, "x_label", "translocation_speed");
    fc_json_array_t *x_range = fc_json_object_set_array(obj, "x_range");
    fc_json_array_append_number(x_range, matrix.x_lower);
    fc_json_array_append_number(x_range, matrix.x_upper);

    fc_json_object_set_string(obj, "y_label", "qscore");
    fc_json_array_t *y_range = fc_json_object_set_array(obj, "y_range");
    fc_json_array_append_number(y_range, 0.0);
    fc_json_array_append_number(y_range, matrix.qscore_bins * matrix.qscore_width);
    fc_json_object_set_value(
        obj, "data", summary_json_heatmap_ndarray(matrix.counts, matrix.qscore_bins, matrix.x_bins));

    fc_json_object_set_value(heatmaps, "translocation_speed_quality", value);
    fc_speed_qscore_2d_matrix_destroy(&matrix);
}

static void summary_json_accuracy_heatmap(
    fc_json_object_t *heatmaps,
    const fc_accuracy_qscore_2d *model,
    size_t barcode)
{
    if (NULL == model) {
        return;
    }
    if (!fc_accuracy_qscore_2d_seen(model, barcode)) {
        fc_json_object_set_null(heatmaps, "alignment_accuracy_quality");
        return;
    }

    fc_accuracy_qscore_2d_matrix matrix = { 0 };
    if (fc_accuracy_qscore_2d_materialize(model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize alignment-accuracy heatmap for summary.json.");
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *obj = fc_json_get_object(value);
    fc_json_object_set_string(obj, "x_label", "alignment_concordance_qscore");
    fc_json_array_t *x_range = fc_json_object_set_array(obj, "x_range");
    fc_json_array_append_number(x_range, matrix.x_lower);
    fc_json_array_append_number(x_range, matrix.x_upper);
    fc_json_object_set_string(obj, "y_label", "qscore");
    fc_json_array_t *y_range = fc_json_object_set_array(obj, "y_range");
    fc_json_array_append_number(y_range, 0.0);
    fc_json_array_append_number(y_range, matrix.qscore_bins * matrix.qscore_width);
    fc_json_object_set_value(
        obj, "data", summary_json_heatmap_ndarray(matrix.counts, matrix.qscore_bins, matrix.x_bins));
    fc_json_object_set_value(heatmaps, "alignment_accuracy_quality", value);
    fc_accuracy_qscore_2d_matrix_destroy(&matrix);
}

static void summary_json_time(
    fc_json_object_t *root,
    const fc_time_bin *model,
    size_t barcode,
    bool include_accuracy)
{
    if (NULL == model || !fc_time_bin_seen(model, barcode)) {
        fc_json_object_set_null(root, "time");
        return;
    }

    fc_time_bin_table table = { 0 };
    if (fc_time_bin_materialize(model, barcode, &table) != 0) {
        FC_LOG_FATAL("Failed to materialize time section for summary.json.");
    }

    if (table.n == 0) {
        fc_time_bin_table_destroy(&table);
        fc_json_object_set_null(root, "time");
        return;
    }

    fc_json_object_t *time = fc_json_object_set_object(root, "time");

    fc_time_bin_skip_counts skipped = fc_time_bin_skipped_counts(model);
    fc_json_object_t *skip = fc_json_object_set_object(time, "skipped");
    fc_json_object_set_number(skip, "missing_start_time", (double)skipped.missing);
    fc_json_object_set_number(skip, "invalid_start_time", (double)skipped.invalid);
    fc_json_object_set_number(skip, "missing_duration", (double)skipped.duration_missing);
    fc_json_object_set_number(skip, "invalid_duration", (double)skipped.duration_invalid);

    // Pack time-series data as columnar arrays
    size_t *read_count = fc_alloc(table.n, sizeof(size_t), "read_count");
    size_t *base_count = fc_alloc(table.n, sizeof(size_t), "base_count");
    size_t *cum_read = fc_alloc(table.n, sizeof(size_t), "cum_read");
    size_t *cum_base = fc_alloc(table.n, sizeof(size_t), "cum_base");
    size_t *accuracy_read_count = fc_alloc(table.n, sizeof(size_t), "accuracy_read_count");
    double *time_start = fc_alloc(table.n, sizeof(double), "time_start");
    double *time_end = fc_alloc(table.n, sizeof(double), "time_end");
    double *len_mean = fc_alloc(table.n, sizeof(double), "len_mean");
    double *len_q05 = fc_alloc(table.n, sizeof(double), "len_q05");
    double *len_q10 = fc_alloc(table.n, sizeof(double), "len_q10");
    double *len_q25 = fc_alloc(table.n, sizeof(double), "len_q25");
    double *len_q50 = fc_alloc(table.n, sizeof(double), "len_q50");
    double *len_q75 = fc_alloc(table.n, sizeof(double), "len_q75");
    double *len_q90 = fc_alloc(table.n, sizeof(double), "len_q90");
    double *len_q95 = fc_alloc(table.n, sizeof(double), "len_q95");
    double *q_mean = fc_alloc(table.n, sizeof(double), "q_mean");
    double *q_q05 = fc_alloc(table.n, sizeof(double), "q_q05");
    double *q_q10 = fc_alloc(table.n, sizeof(double), "q_q10");
    double *q_q25 = fc_alloc(table.n, sizeof(double), "q_q25");
    double *q_q50 = fc_alloc(table.n, sizeof(double), "q_q50");
    double *q_q75 = fc_alloc(table.n, sizeof(double), "q_q75");
    double *q_q90 = fc_alloc(table.n, sizeof(double), "q_q90");
    double *q_q95 = fc_alloc(table.n, sizeof(double), "q_q95");
    double *spd_mean = fc_alloc(table.n, sizeof(double), "spd_mean");
    double *spd_q05 = fc_alloc(table.n, sizeof(double), "spd_q05");
    double *spd_q10 = fc_alloc(table.n, sizeof(double), "spd_q10");
    double *spd_q25 = fc_alloc(table.n, sizeof(double), "spd_q25");
    double *spd_q50 = fc_alloc(table.n, sizeof(double), "spd_q50");
    double *spd_q75 = fc_alloc(table.n, sizeof(double), "spd_q75");
    double *spd_q90 = fc_alloc(table.n, sizeof(double), "spd_q90");
    double *spd_q95 = fc_alloc(table.n, sizeof(double), "spd_q95");
    double *accuracy_mean = fc_alloc(table.n, sizeof(double), "accuracy_mean");
    double *accuracy_q05 = fc_alloc(table.n, sizeof(double), "accuracy_q05");
    double *accuracy_q10 = fc_alloc(table.n, sizeof(double), "accuracy_q10");
    double *accuracy_q25 = fc_alloc(table.n, sizeof(double), "accuracy_q25");
    double *accuracy_q50 = fc_alloc(table.n, sizeof(double), "accuracy_q50");
    double *accuracy_q75 = fc_alloc(table.n, sizeof(double), "accuracy_q75");
    double *accuracy_q90 = fc_alloc(table.n, sizeof(double), "accuracy_q90");
    double *accuracy_q95 = fc_alloc(table.n, sizeof(double), "accuracy_q95");

    for (size_t i = 0; i < table.n; i++) {
        const fc_time_bin_row *row = &table.rows[i];
        time_start[i] = row->time_bin_start_s;
        time_end[i] = row->time_bin_end_s;
        read_count[i] = row->read_count;
        base_count[i] = row->base_count;
        cum_read[i] = row->cumulative_read_count;
        cum_base[i] = row->cumulative_base_count;
        len_mean[i] = row->length_mean;
        len_q05[i] = row->length_quantiles.q05;
        len_q10[i] = row->length_quantiles.q10;
        len_q25[i] = row->length_quantiles.q25;
        len_q50[i] = row->length_quantiles.q50;
        len_q75[i] = row->length_quantiles.q75;
        len_q90[i] = row->length_quantiles.q90;
        len_q95[i] = row->length_quantiles.q95;
        q_mean[i] = row->qscore_mean;
        q_q05[i] = row->qscore_quantiles.q05;
        q_q10[i] = row->qscore_quantiles.q10;
        q_q25[i] = row->qscore_quantiles.q25;
        q_q50[i] = row->qscore_quantiles.q50;
        q_q75[i] = row->qscore_quantiles.q75;
        q_q90[i] = row->qscore_quantiles.q90;
        q_q95[i] = row->qscore_quantiles.q95;
        spd_mean[i] = row->speed_mean;
        spd_q05[i] = row->speed_quantiles.q05;
        spd_q10[i] = row->speed_quantiles.q10;
        spd_q25[i] = row->speed_quantiles.q25;
        spd_q50[i] = row->speed_quantiles.q50;
        spd_q75[i] = row->speed_quantiles.q75;
        spd_q90[i] = row->speed_quantiles.q90;
        spd_q95[i] = row->speed_quantiles.q95;
        accuracy_read_count[i] = row->accuracy_read_count;
        accuracy_mean[i] = row->accuracy_mean;
        accuracy_q05[i] = row->accuracy_quantiles.q05;
        accuracy_q10[i] = row->accuracy_quantiles.q10;
        accuracy_q25[i] = row->accuracy_quantiles.q25;
        accuracy_q50[i] = row->accuracy_quantiles.q50;
        accuracy_q75[i] = row->accuracy_quantiles.q75;
        accuracy_q90[i] = row->accuracy_quantiles.q90;
        accuracy_q95[i] = row->accuracy_quantiles.q95;
    }

    // Keep timing skip counts separate from the self-describing table payload.
    fc_json_object_t *data = fc_json_object_set_object(time, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(
        data, table.n, include_accuracy ? 39 : 30);
    fc_json_array_append_value(columns, fc_json_column_f64("time_bin_start_s", time_start, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("time_bin_end_s", time_end, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("read_count", read_count, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("base_count", base_count, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("cumulative_read_count", cum_read, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("cumulative_base_count", cum_base, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_mean", len_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q05", len_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q10", len_q10, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q25", len_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q50", len_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q75", len_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q90", len_q90, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q95", len_q95, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_mean", q_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q05", q_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q10", q_q10, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q25", q_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q50", q_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q75", q_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q90", q_q90, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q95", q_q95, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_mean", spd_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q05", spd_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q10", spd_q10, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q25", spd_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q50", spd_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q75", spd_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q90", spd_q90, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q95", spd_q95, table.n));
    if (include_accuracy) {
        fc_json_array_append_value(columns, fc_json_column_uint_adaptive("accuracy_read_count", accuracy_read_count, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_mean", accuracy_mean, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q05", accuracy_q05, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q10", accuracy_q10, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q25", accuracy_q25, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q50", accuracy_q50, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q75", accuracy_q75, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q90", accuracy_q90, table.n));
        fc_json_array_append_value(columns, fc_json_column_f64("accuracy_q95", accuracy_q95, table.n));
    }

    free(time_start);
    free(time_end);
    free(read_count);
    free(base_count);
    free(cum_read);
    free(cum_base);
    free(accuracy_read_count);
    free(len_mean);
    free(len_q05);
    free(len_q10);
    free(len_q25);
    free(len_q50);
    free(len_q75);
    free(len_q90);
    free(len_q95);
    free(q_mean);
    free(q_q05);
    free(q_q10);
    free(q_q25);
    free(q_q50);
    free(q_q75);
    free(q_q90);
    free(q_q95);
    free(spd_mean);
    free(spd_q05);
    free(spd_q10);
    free(spd_q25);
    free(spd_q50);
    free(spd_q75);
    free(spd_q90);
    free(spd_q95);
    free(accuracy_mean);
    free(accuracy_q05);
    free(accuracy_q10);
    free(accuracy_q25);
    free(accuracy_q50);
    free(accuracy_q75);
    free(accuracy_q90);
    free(accuracy_q95);

    fc_time_bin_table_destroy(&table);
}

static void summary_json_channels(fc_json_object_t *root, const fc_channel_summary *model, size_t barcode)
{
    if (NULL == model || !fc_channel_summary_seen(model, barcode)) {
        fc_json_object_set_null(root, "channels");
        return;
    }

    fc_channel_summary_table table = { 0 };
    if (fc_channel_summary_materialize_sparse(model, barcode, &table) != 0) {
        FC_LOG_FATAL("Failed to materialize channel section for summary.json.");
    }

    if (table.n == 0) {
        fc_channel_summary_table_destroy(&table);
        fc_json_object_set_null(root, "channels");
        return;
    }

    fc_json_object_t *channels = fc_json_object_set_object(root, "channels");
    fc_json_object_set_string(channels, "device", summary_channel_device_name(table.highest_channel));

    fc_channel_summary_skip_counts skipped = fc_channel_summary_skipped_counts(model);
    fc_json_object_t *skip = fc_json_object_set_object(channels, "skipped");
    fc_json_object_set_number(skip, "missing_channel", (double)skipped.missing_channel);
    fc_json_object_set_number(skip, "invalid_channel", (double)skipped.invalid_channel);
    fc_json_object_set_number(skip, "missing_mux", (double)skipped.missing_mux);
    fc_json_object_set_number(skip, "invalid_mux", (double)skipped.invalid_mux);

    size_t *channel_arr = fc_alloc(table.n, sizeof(size_t), "channel array");
    size_t *mux_arr = fc_alloc(table.n, sizeof(size_t), "mux array");
    size_t *read_count_arr = fc_alloc(table.n, sizeof(size_t), "read_count array");
    size_t *base_count_arr = fc_alloc(table.n, sizeof(size_t), "base_count array");
    size_t *quality_read_count_arr = fc_alloc(table.n, sizeof(size_t), "quality_read_count array");
    double *mean_quality_arr = fc_alloc(table.n, sizeof(double), "mean_quality array");
    size_t *concordance_read_count_arr = fc_alloc(
        table.n, sizeof(size_t), "concordance_read_count array");
    double *mean_concordance_qscore_arr = fc_alloc(
        table.n, sizeof(double), "mean_concordance_qscore array");

    for (size_t i = 0; i < table.n; i++) {
        const fc_channel_summary_row *row = &table.rows[i];
        channel_arr[i] = row->channel;
        mux_arr[i] = row->mux;
        read_count_arr[i] = row->read_count;
        base_count_arr[i] = row->base_count;
        quality_read_count_arr[i] = row->quality_read_count;
        mean_quality_arr[i] = row->mean_quality;
        concordance_read_count_arr[i] = row->concordance_read_count;
        mean_concordance_qscore_arr[i] = row->mean_concordance_qscore;
    }

    fc_json_object_t *data = fc_json_object_set_object(channels, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(data, table.n, 8);
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("channel", channel_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("mux", mux_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("read_count", read_count_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("base_count", base_count_arr, table.n));
    fc_json_array_append_value(
        columns, fc_json_column_uint_adaptive("quality_read_count", quality_read_count_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("mean_quality", mean_quality_arr, table.n));
    fc_json_array_append_value(
        columns,
        fc_json_column_uint_adaptive("concordance_read_count", concordance_read_count_arr, table.n));
    fc_json_array_append_value(
        columns, fc_json_column_f64("mean_concordance_qscore", mean_concordance_qscore_arr, table.n));

    free(channel_arr);
    free(mux_arr);
    free(read_count_arr);
    free(base_count_arr);
    free(quality_read_count_arr);
    free(mean_quality_arr);
    free(concordance_read_count_arr);
    free(mean_concordance_qscore_arr);

    fc_channel_summary_table_destroy(&table);
}

static void summary_json_flagstats(
    fc_json_object_t *alignment,
    const summary_json_sink *sink,
    size_t barcode)
{
    if (!sink->cfg.bam_input) {
        return;
    }

    const fc_alignment_counts *counts = fc_alignment_summary_for_output(
        sink->cfg.alignment_summary,
        barcode);
    if (NULL == counts || !counts->seen) {
        fc_json_object_t *record_counts = fc_json_object_set_object(alignment, "record_counts");
        size_t zero_counts[FC_FLAGSTAT_COUNT] = { 0 };
        for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
            fc_json_object_set_number(record_counts, fc_flagstat_column_names[col], (double)zero_counts[col]);
        }
        fc_json_object_set_number(alignment, "mapped_record_count", 0.0);
        fc_json_object_set_number(alignment, "mapped_record_percent", 0.0);
        fc_json_object_set_number(alignment, "primary_mapped_record_count", 0.0);
        fc_json_object_set_number(alignment, "primary_mapped_record_percent", 0.0);
        fc_json_object_set_number(alignment, "read_count", 0.0);
        fc_json_object_set_number(alignment, "mapped_read_count", 0.0);
        fc_json_object_set_number(alignment, "mapped_read_percent", 0.0);
        fc_json_object_set_array(alignment, "per_reference_counts");
        return;
    }

    fc_alignment_totals totals = { 0 };
    fc_alignment_counts_materialize(counts, &totals);

    fc_json_object_t *record_counts = fc_json_object_set_object(alignment, "record_counts");
    for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
        fc_json_object_set_number(
            record_counts,
            fc_flagstat_column_names[col],
            (double)totals.record_counts[col]);
    }

    fc_json_object_set_number(
        alignment,
        "mapped_record_count",
        (double)totals.mapped_record_count);
    fc_json_object_set_number(
        alignment,
        "mapped_record_percent",
        totals.record_counts[FC_FLAGSTAT_TOTAL] > 0
            ? 100.0 * totals.mapped_record_count / totals.record_counts[FC_FLAGSTAT_TOTAL]
            : 0.0);
    fc_json_object_set_number(
        alignment,
        "primary_mapped_record_count",
        (double)totals.primary_mapped_record_count);
    fc_json_object_set_number(
        alignment,
        "primary_mapped_record_percent",
        totals.record_counts[FC_FLAGSTAT_TOTAL] > 0
            ? 100.0 * totals.primary_mapped_record_count / totals.record_counts[FC_FLAGSTAT_TOTAL]
            : 0.0);
    fc_json_object_set_number(alignment, "read_count", (double)totals.read_count);
    fc_json_object_set_number(alignment, "mapped_read_count", (double)totals.mapped_read_count);
    fc_json_object_set_number(alignment, "mapped_read_percent", 100.0 * totals.mapping_rate);

    fc_json_array_t *rows = fc_json_object_set_array(alignment, "per_reference_counts");
    if (sink->cfg.region_mode) {
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *obj = fc_json_get_object(value);
        fc_json_object_set_string(obj, "ref", sink->cfg.region_label ? sink->cfg.region_label : "region");
        const size_t *row_counts = counts->reference_counts[0];
        for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
            fc_json_object_set_number(
                obj,
                fc_flagstat_column_names[col],
                (double)row_counts[col]);
        }
        fc_json_array_append_value(rows, value);
    } else {
        for (size_t i = 0; i < counts->n_refs; i++) {
            fc_json_value_t *value = fc_json_object_value();
            fc_json_object_t *obj = fc_json_get_object(value);
            const char *label = sink->cfg.hdr && i < (size_t)sink->cfg.hdr->n_targets
                ? sam_hdr_tid2name((sam_hdr_t *)sink->cfg.hdr, (int)i)
                : NULL;
            char tmp[32] = { 0 };
            if (NULL == label) {
                FC_SNPRINTF_CHECK(tmp, sizeof(tmp), "%zu", i);
                label = tmp;
            }
            const size_t *row_counts = counts->reference_counts[i];
            fc_json_object_set_string(obj, "ref", label);
            for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
                fc_json_object_set_number(
                    obj,
                    fc_flagstat_column_names[col],
                    (double)row_counts[col]);
            }
            fc_json_array_append_value(rows, value);
        }
    }
    if (NULL != counts->unmapped_counts) {
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *obj = fc_json_get_object(value);
        fc_json_object_set_string(obj, "ref", "*");
        for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
            fc_json_object_set_number(
                obj,
                fc_flagstat_column_names[col],
                (double)counts->unmapped_counts[col]);
        }
        fc_json_array_append_value(rows, value);
    }
}

static void summary_json_alignment_distributions(
    fc_json_object_t *alignment,
    const fc_stats_engine *engine)
{
    fc_json_object_t *distributions = fc_json_object_set_object(alignment, "distributions");

    fc_json_value_t *value = summary_json_fixed_distribution_value(engine ? engine->acc_stats : NULL);
    fc_json_object_set_value(distributions, "accuracy", value);

    value = summary_json_fixed_distribution_value(engine ? engine->cov_stats : NULL);
    fc_json_object_set_value(distributions, "query_coverage", value);

    value = summary_json_length_distribution_value(
        engine ? engine->polya_stats : NULL, 0, 0, false, false);
    fc_json_object_set_value(distributions, "poly_a_length", value);
}

static void summary_json_alignment(
    fc_json_object_t *root,
    const summary_json_sink *sink,
    size_t barcode,
    const fc_stats_engine *engine)
{
    if (!sink->cfg.bam_input) {
        return;
    }

    fc_json_object_t *alignment = fc_json_object_set_object(root, "alignment");
    summary_json_flagstats(alignment, sink, barcode);
    summary_json_alignment_distributions(alignment, engine);
}

static void summary_json_coverage_table(
    fc_json_object_t *object, const fc_coverage_summary_target *target)
{
    size_t n_rows = target->n_summary_rows;
    const char **chrom = n_rows ? fc_alloc(n_rows, sizeof(*chrom), "coverage JSON summary chrom") : NULL;
    size_t *start = n_rows ? fc_alloc(n_rows, sizeof(*start), "coverage JSON summary start") : NULL;
    size_t *end = n_rows ? fc_alloc(n_rows, sizeof(*end), "coverage JSON summary end") : NULL;
    size_t *length = n_rows ? fc_alloc(n_rows, sizeof(*length), "coverage JSON summary length") : NULL;
    size_t *bases = n_rows ? fc_alloc(n_rows, sizeof(*bases), "coverage JSON summary bases") : NULL;
    double *mean = n_rows ? fc_alloc(n_rows, sizeof(*mean), "coverage JSON summary mean") : NULL;
    size_t *min = n_rows ? fc_alloc(n_rows, sizeof(*min), "coverage JSON summary min") : NULL;
    size_t *max = n_rows ? fc_alloc(n_rows, sizeof(*max), "coverage JSON summary max") : NULL;

    for (size_t row = 0; row < n_rows; ++row) {
        const fc_coverage_summary_row *source = &target->summary_rows[row];
        chrom[row] = source->chrom;
        start[row] = source->start;
        end[row] = source->end;
        length[row] = source->length;
        bases[row] = source->bases;
        mean[row] = source->mean;
        min[row] = source->min;
        max[row] = source->max;
    }

    fc_json_object_t *data = fc_json_object_set_object(object, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(
        data, n_rows, 8 + target->n_thresholds);
    fc_json_array_append_value(columns, fc_json_column_string("chrom", chrom, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("start", start, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("end", end, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("length", length, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("bases", bases, n_rows));
    fc_json_array_append_value(columns, fc_json_column_f64("mean", mean, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("min", min, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("max", max, n_rows));
    for (size_t threshold = 0; threshold < target->n_thresholds; ++threshold) {
        char name[32];
        FC_SNPRINTF_CHECK(name, sizeof(name), "%ux", target->thresholds[threshold]);
        double *fractions = n_rows
            ? fc_alloc(n_rows, sizeof(*fractions), "coverage JSON threshold fractions")
            : NULL;
        for (size_t row = 0; row < n_rows; ++row) {
            fractions[row] = target->summary_rows[row].threshold_fractions[threshold];
        }
        fc_json_array_append_value(columns, fc_json_column_f64(name, fractions, n_rows));
        free(fractions);
    }

    free(chrom);
    free(start);
    free(end);
    free(length);
    free(bases);
    free(mean);
    free(min);
    free(max);
}

static void summary_json_coverage_distribution(
    fc_json_object_t *object, const fc_coverage_summary_target *target)
{
    size_t n_rows = target->n_distribution_rows;
    const char **chrom = n_rows
        ? fc_alloc(n_rows, sizeof(*chrom), "coverage JSON distribution chrom")
        : NULL;
    size_t *depth = n_rows ? fc_alloc(n_rows, sizeof(*depth), "coverage JSON distribution depth") : NULL;
    double *fraction = n_rows
        ? fc_alloc(n_rows, sizeof(*fraction), "coverage JSON distribution fraction")
        : NULL;
    for (size_t row = 0; row < n_rows; ++row) {
        chrom[row] = target->distribution_rows[row].chrom;
        depth[row] = target->distribution_rows[row].depth;
        fraction[row] = target->distribution_rows[row].fraction;
    }

    fc_json_object_t *data = fc_json_object_set_object(object, "data");
    fc_json_array_t *columns = fc_json_object_set_dataframe(data, n_rows, 3);
    fc_json_array_append_value(columns, fc_json_column_string("chrom", chrom, n_rows));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("depth", depth, n_rows));
    fc_json_array_append_value(columns, fc_json_column_f64("fraction", fraction, n_rows));

    free(chrom);
    free(depth);
    free(fraction);
}

static void summary_json_coverage(
    fc_json_object_t *root, const fc_coverage_summaries *summaries, size_t barcode)
{
    const fc_coverage_summary *summary = fc_coverage_summaries_get(summaries, barcode);
    if (NULL == summary) {
        return;
    }
    fc_json_object_t *coverage = fc_json_object_set_object(root, "coverage");
    fc_json_array_t *targets = fc_json_object_set_array(coverage, "targets");
    for (size_t i = 0; i < summary->n_targets; ++i) {
        const fc_coverage_summary_target *target = &summary->targets[i];
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *target_object = fc_json_get_object(value);
        fc_json_object_set_string(target_object, "name", target->name);
        fc_json_array_t *thresholds = fc_json_object_set_array(target_object, "thresholds");
        for (size_t threshold = 0; threshold < target->n_thresholds; ++threshold) {
            fc_json_array_append_number(thresholds, (double)target->thresholds[threshold]);
        }
        fc_json_object_t *table = fc_json_object_set_object(target_object, "summary");
        summary_json_coverage_table(table, target);
        fc_json_object_t *distribution = fc_json_object_set_object(target_object, "distribution");
        summary_json_coverage_distribution(distribution, target);
        fc_json_array_append_value(targets, value);
    }
}

static FILE *summary_json_open_output(const summary_json_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (NULL == dir) {
        FC_LOG_FATAL("Failed to create output directory for summary.json (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "summary.json");
    free(dir);
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (NULL == fp) {
        FC_LOG_FATAL("Cannot open summary.json output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

static void summary_json_write_barcode_output(summary_json_sink *sink, size_t barcode)
{
    summary_json_barcode *state = &sink->barcodes[barcode];
    size_t engine_idx = sink->cfg.use_barcodes ? barcode : 0;
    const fc_stats_engine *engine = engine_idx < sink->cfg.engine_count ? sink->cfg.engines[engine_idx] : NULL;
    const read_stats *length_stats = summary_read_length_stats(engine);
    const read_stats *quality_stats = summary_quality_stats(engine);

    char barcode_buf[32] = { 0 };
    const char *barcode_id = sink->cfg.use_barcodes
        ? fc_sink_barcode_stem(barcode, barcode_buf, sizeof(barcode_buf))
        : "all";

    fc_run_summary_row row = { 0 };
    if (fc_run_summary_materialize(
            sink->cfg.run_summary,
            barcode,
            engine,
            sink->cfg.bam_input,
            length_stats,
            sink->cfg.sample,
            barcode_id,
            &row)
        != 0) {
        FC_LOG_FATAL("Failed to materialize run_summary for summary.json barcode %zu.", barcode);
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *root = fc_json_get_object(value);
    fc_json_object_set_string(root, "schema_id", FC_SUMMARY_JSON_SCHEMA_ID);
    fc_json_object_set_number(root, "schema_version", FC_SUMMARY_JSON_SCHEMA_VERSION);

    fc_json_object_t *summary = fc_json_object_set_object(root, "run_summary");
    summary_json_set_string_or_null(summary, "sample_id", row.sample_id);
    summary_json_set_string_or_null(summary, "input_type", row.input_type);
    summary_json_set_string_or_null(summary, "barcode_id", row.barcode_id);
    summary_json_set_string_or_null(summary, "barcode_alias", row.barcode_alias);
    summary_json_set_semicolon_list_or_null(summary, "run_ids", row.run_id);
    summary_json_set_semicolon_list_or_null(summary, "basecaller_models", row.basecaller_model);
    summary_json_set_semicolon_list_or_null(summary, "modbase_models", row.modbase_models);
    summary_json_set_semicolon_list_or_null(summary, "flow_cell_ids", row.flow_cell_id);
    summary_json_set_string_or_null(summary, "run_date", row.run_date);
    fc_json_object_set_number(summary, "run_duration_s", row.run_duration_s);
    fc_json_object_set_number(summary, "read_count", (double)row.read_count);
    fc_json_object_set_number(summary, "base_count", (double)row.base_count);
    fc_json_object_set_number(summary, "records_seen", (double)row.records_seen);
    fc_json_object_set_number(summary, "records_passed_filters", (double)row.records_passed_filters);
    fc_json_object_set_number(summary, "records_too_long", (double)row.records_too_long);
    fc_json_object_set_number(summary, "records_too_short", (double)row.records_too_short);
    fc_json_object_set_number(summary, "records_low_quality", (double)row.records_low_quality);
    fc_json_object_set_number(summary, "records_dust_masked", (double)row.records_dust_masked);
    fc_json_object_set_number(summary, "n5_bp", (double)row.n5_bp);
    fc_json_object_set_number(summary, "n10_bp", (double)row.n10_bp);
    fc_json_object_set_number(summary, "n25_bp", (double)row.n25_bp);
    fc_json_object_set_number(summary, "n50_bp", (double)row.n50_bp);
    fc_json_object_set_number(summary, "n75_bp", (double)row.n75_bp);
    fc_json_object_set_number(summary, "n90_bp", (double)row.n90_bp);
    fc_json_object_set_number(summary, "n95_bp", (double)row.n95_bp);
    fc_json_object_set_number(summary, "l5", (double)row.l5);
    fc_json_object_set_number(summary, "l10", (double)row.l10);
    fc_json_object_set_number(summary, "l25", (double)row.l25);
    fc_json_object_set_number(summary, "l50", (double)row.l50);
    fc_json_object_set_number(summary, "l75", (double)row.l75);
    fc_json_object_set_number(summary, "l90", (double)row.l90);
    fc_json_object_set_number(summary, "l95", (double)row.l95);

    fc_json_object_t *sample = fc_json_object_set_object(root, "sample");
    summary_json_set_string_or_null(sample, "sample_id", row.sample_id);
    summary_json_set_string_or_null(sample, "input_type", row.input_type);
    summary_json_set_string_or_null(sample, "barcode_id", row.barcode_id);
    summary_json_set_string_or_null(sample, "barcode_alias", row.barcode_alias);

    summary_json_files(root, state);
    summary_json_counts(root, &row, state);
    summary_json_channels(root, sink->cfg.channel_summary, barcode);

    fc_json_object_t *distributions = fc_json_object_set_object(root, "distributions");
    fc_json_value_t *length_distribution = summary_json_length_distribution_value(
        length_stats, row.base_count, row.read_count, true, true);
    fc_json_object_set_value(distributions, "read_length", length_distribution);

    fc_json_value_t *quality_distribution = summary_json_fixed_distribution_value(quality_stats);
    fc_json_object_set_value(distributions, "quality", quality_distribution);

    summary_json_time(root, sink->cfg.time_bin, barcode, sink->cfg.accuracy_qscore != NULL);

    fc_json_object_t *heatmaps = fc_json_object_set_object(root, "heatmaps");
    summary_json_length_heatmap(heatmaps, sink->cfg.length_qscore, barcode);
    summary_json_speed_heatmap(heatmaps, sink->cfg.speed_qscore, barcode);
    summary_json_accuracy_heatmap(heatmaps, sink->cfg.accuracy_qscore, barcode);
    summary_json_alignment(root, sink, barcode, engine);
    summary_json_coverage(root, sink->cfg.coverage_summaries, barcode);

    FILE *fp = summary_json_open_output(sink, barcode);
    char *json = fc_json_serialize_to_string(value);
    fputs(json, fp);
    fputc('\n', fp);
    fclose(fp);
    free(json);
    fc_json_value_free(value);
    fc_run_summary_row_destroy(&row);
}

static void summary_json_sink_begin_file(void *ctx, const char *fname)
{
    summary_json_sink *sink = (summary_json_sink *)ctx;
    (void)fname;
    for (size_t i = 0; i < sink->cfg.barcode_count; i++) {
        summary_file_accum_reset(&sink->barcodes[i].current);
    }
}

static void summary_json_sink_add_read(
    void *ctx,
    const fc_read *read,
    const fc_read_ctx *read_ctx,
    const char *fname,
    const sam_hdr_t *hdr)
{
    summary_json_sink *sink = (summary_json_sink *)ctx;
    (void)fname;
    (void)hdr;
    if (!fc_is_canonical_read_record(read->bam)) {
        return;
    }

    size_t barcode = (sink->cfg.use_barcodes && read_ctx) ? read_ctx->barcode : 0;
    if (barcode >= sink->cfg.barcode_count) {
        return;
    }
    summary_json_barcode *state = &sink->barcodes[barcode];
    summary_barcode_add_read(state, read);
}

static void summary_json_sink_end_file(void *ctx, const char *fname, const fc_file_stats *stats)
{
    summary_json_sink *sink = (summary_json_sink *)ctx;
    if (!sink->cfg.use_barcodes) {
        summary_barcode_append_file(&sink->barcodes[0], fname, stats, true);
        return;
    }
    for (size_t i = 0; i < sink->cfg.barcode_count; i++) {
        summary_barcode_append_file(&sink->barcodes[i], fname, NULL, false);
    }
}

static void summary_json_sink_finish(void *ctx)
{
    summary_json_sink *sink = (summary_json_sink *)ctx;
    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (sink->cfg.run_summary && fc_run_summary_seen(sink->cfg.run_summary, b)) {
                summary_json_write_barcode_output(sink, b);
            }
        }
    } else {
        summary_json_write_barcode_output(sink, 0);
    }
}

static void summary_json_sink_destroy(void *ctx)
{
    summary_json_sink *sink = (summary_json_sink *)ctx;
    if (NULL == sink) {
        return;
    }
    for (size_t i = 0; i < sink->cfg.barcode_count; i++) {
        summary_barcode_destroy(&sink->barcodes[i]);
    }
    free(sink->barcodes);
    free(sink->output_root);
    free(sink->sample);
    free(sink);
}

fc_sink fc_summary_json_sink_create(const fc_summary_json_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (NULL == opts
        || NULL == opts->output_root
        || NULL == opts->run_summary
        || NULL == opts->channel_summary
        || NULL == opts->length_qscore
        || NULL == opts->speed_qscore
        || NULL == opts->time_bin) {
        return sink;
    }

    summary_json_sink *state = fc_alloc(1, sizeof(*state), "summary_json sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);
    if (NULL != opts->sample) {
        state->sample = fc_strdup(opts->sample, "summary_json sample");
        state->cfg.sample = state->sample;
    }
    state->barcodes = fc_alloc(state->cfg.barcode_count, sizeof(*state->barcodes), "summary_json barcodes");

    sink.ctx = state;
    sink.begin_file = summary_json_sink_begin_file;
    sink.add_read = summary_json_sink_add_read;
    sink.end_file = summary_json_sink_end_file;
    sink.finish = summary_json_sink_finish;
    sink.destroy = summary_json_sink_destroy;
    return sink;
}
