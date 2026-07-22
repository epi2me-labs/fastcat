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
#define FC_SUMMARY_JSON_FLAG_COLS 9

typedef struct {
    size_t n_refs;
    size_t **counts;
    size_t *unmapped;
} summary_flag_stats;

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
    summary_flag_stats *flagstats;
} summary_json_barcode;

typedef struct {
    fc_summary_json_sink_opts cfg;
    char *output_root;
    char *sample;
    size_t n_refs;
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

static void summary_flag_stats_destroy(summary_flag_stats *stats)
{
    if (NULL == stats) {
        return;
    }
    for (size_t i = 0; i < stats->n_refs; i++) {
        free(stats->counts[i]);
    }
    free(stats->counts);
    free(stats->unmapped);
    free(stats);
}

static summary_flag_stats *summary_flag_stats_create(size_t n_refs, bool store_unmapped)
{
    summary_flag_stats *stats = fc_alloc(1, sizeof(*stats), "summary_json flagstats");
    stats->n_refs = n_refs;
    stats->counts = n_refs > 0 ? fc_alloc(n_refs, sizeof(*stats->counts), "summary_json flagstats counts") : NULL;
    stats->unmapped = store_unmapped
        ? fc_alloc(FC_SUMMARY_JSON_FLAG_COLS, sizeof(*stats->unmapped), "summary_json flagstats unmapped")
        : NULL;
    for (size_t i = 0; i < n_refs; i++) {
        stats->counts[i] = fc_alloc(
            FC_SUMMARY_JSON_FLAG_COLS,
            sizeof(*stats->counts[i]),
            "summary_json flagstats row");
    }
    return stats;
}

static void summary_flagstats_add_counts(const bam1_t *bam, size_t *counts, int duplex_code)
{
    counts[0] += 1;
    counts[1] += fc_is_primary(bam);
    counts[2] += ((bam->core.flag & BAM_FSECONDARY) != 0);
    counts[3] += ((bam->core.flag & BAM_FSUPPLEMENTARY) != 0);
    counts[4] += ((bam->core.flag & BAM_FUNMAP) != 0);
    counts[5] += ((bam->core.flag & BAM_FQCFAIL) != 0);
    counts[6] += ((bam->core.flag & BAM_FDUP) != 0);
    counts[7] += (duplex_code == 1);
    counts[8] += (duplex_code == -1);
}

static summary_flag_stats *summary_ensure_flagstats(summary_json_sink *sink, size_t barcode)
{
    if (barcode >= sink->cfg.barcode_count) {
        return NULL;
    }
    if (NULL == sink->barcodes[barcode].flagstats) {
        sink->barcodes[barcode].flagstats = summary_flag_stats_create(sink->n_refs, true);
    }
    return sink->barcodes[barcode].flagstats;
}

static void summary_counter_add(kh_counter_t **counter, const char *value)
{
    if (NULL == value || '\0' == value[0]) {
        return;
    }
    if (NULL == *counter) {
        *counter = kh_counter_init();
    }
    kh_counter_increment(*counter, (char *)value);
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
    summary_flag_stats_destroy(state->flagstats);
    memset(state, 0, sizeof(*state));
}

static const read_stats *summary_read_length_stats(const fc_stats_engine *engine, bool unmapped_override)
{
    if (NULL == engine) {
        return NULL;
    }
    if (unmapped_override) {
        return engine->length_stats_unmapped;
    }
    double aligned = 0.0;
    double unmapped = 0.0;
    if (NULL != engine->length_stats) {
        for (size_t i = 0; i < engine->length_stats->n; i++) {
            aligned += (double)engine->length_stats->counts[i];
        }
    }
    if (NULL != engine->length_stats_unmapped) {
        for (size_t i = 0; i < engine->length_stats_unmapped->n; i++) {
            unmapped += (double)engine->length_stats_unmapped->counts[i];
        }
    }
    if (aligned > 0.0) {
        return engine->length_stats;
    }
    if (unmapped > 0.0) {
        return engine->length_stats_unmapped;
    }
    return engine->length_stats;
}

static const read_stats *summary_quality_stats(
    const fc_stats_engine *engine,
    bool bam_input,
    bool unmapped_override)
{
    if (NULL == engine) {
        return NULL;
    }
    if (unmapped_override) {
        return engine->qual_stats_unmapped;
    }

    double aligned = 0.0;
    double unmapped = 0.0;
    if (NULL != engine->qual_stats) {
        for (size_t i = 0; i < engine->qual_stats->n; i++) {
            aligned += (double)engine->qual_stats->counts[i];
        }
    }
    if (NULL != engine->qual_stats_unmapped) {
        for (size_t i = 0; i < engine->qual_stats_unmapped->n; i++) {
            unmapped += (double)engine->qual_stats_unmapped->counts[i];
        }
    }
    if (aligned > 0.0) {
        return engine->qual_stats;
    }
    if (bam_input && unmapped > 0.0) {
        return engine->qual_stats_unmapped;
    }
    return engine->qual_stats;
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

static void summary_json_length_distribution(fc_json_object_t *distributions, const read_stats *stats)
{
    if (NULL == stats) {
        fc_json_object_set_null(distributions, "read_length");
        return;
    }

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *hist = fc_json_get_object(value);

    // Collect non-zero bins
    size_t *starts_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length starts");
    size_t *ends_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length ends");
    size_t *counts_arr = fc_alloc(stats->n, sizeof(size_t), "summary_json length counts");
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

    // Encode as columnar table (consistent with other data tables)
    fc_json_object_t *data = fc_json_object_set_object(hist, "__data__");
    fc_json_array_t *shape = fc_json_object_set_array(data, "shape");
    fc_json_array_append_number(shape, (double)n_bins);

    fc_json_array_t *columns = fc_json_object_set_array(data, "columns");
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("start", starts_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("end", ends_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("count", counts_arr, n_bins));

    free(starts_arr);
    free(ends_arr);
    free(counts_arr);

    fc_json_object_set_number(hist, "total_count", total);
    if (overflow_count > 0) {
        fc_json_object_set_number(hist, "overflow_start_bp", (double)overflow_start);
        fc_json_object_set_number(hist, "overflow_count", (double)overflow_count);
    } else {
        fc_json_object_set_null(hist, "overflow_start_bp");
        fc_json_object_set_null(hist, "overflow_count");
    }
    fc_json_object_set_value(distributions, "read_length", value);
}

static void summary_json_quality_distribution(fc_json_object_t *distributions, const read_stats *stats)
{
    if (NULL == stats) {
        fc_json_object_set_null(distributions, "quality");
        return;
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

    // Encode as columnar table (consistent with other data tables)
    fc_json_object_t *data = fc_json_object_set_object(hist, "__data__");
    fc_json_array_t *shape = fc_json_object_set_array(data, "shape");
    fc_json_array_append_number(shape, (double)n_bins);

    fc_json_array_t *columns = fc_json_object_set_array(data, "columns");
    fc_json_array_append_value(columns, fc_json_column_f64("start", starts_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_f64("end", ends_arr, n_bins));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("count", counts_arr, n_bins));

    free(starts_arr);
    free(ends_arr);
    free(counts_arr);

    fc_json_object_set_value(distributions, "quality", value);
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
        obj,
        "__data__",
        summary_json_heatmap_ndarray(matrix.counts, matrix.qscore_bins, matrix.x_bins));
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
        obj,
        "__data__",
        summary_json_heatmap_ndarray(matrix.counts, matrix.qscore_bins, matrix.x_bins));
    fc_json_object_set_value(heatmaps, "translocation_speed_quality", value);
    fc_speed_qscore_2d_matrix_destroy(&matrix);
}

static void summary_json_time(fc_json_object_t *root, const fc_time_bin *model, size_t barcode)
{
    if (NULL == model || !fc_time_bin_seen(model, barcode)) {
        fc_json_object_set_null(root, "time");
        return;
    }

    fc_time_bin_table table = { 0 };
    if (fc_time_bin_materialize(model, barcode, &table) != 0) {
        FC_LOG_FATAL("Failed to materialize time section for summary.json.");
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
    double *time_start = fc_alloc(table.n, sizeof(double), "time_start");
    double *time_end = fc_alloc(table.n, sizeof(double), "time_end");
    double *len_mean = fc_alloc(table.n, sizeof(double), "len_mean");
    double *len_q05 = fc_alloc(table.n, sizeof(double), "len_q05");
    double *len_q25 = fc_alloc(table.n, sizeof(double), "len_q25");
    double *len_q50 = fc_alloc(table.n, sizeof(double), "len_q50");
    double *len_q75 = fc_alloc(table.n, sizeof(double), "len_q75");
    double *len_q95 = fc_alloc(table.n, sizeof(double), "len_q95");
    double *q_mean = fc_alloc(table.n, sizeof(double), "q_mean");
    double *q_q05 = fc_alloc(table.n, sizeof(double), "q_q05");
    double *q_q25 = fc_alloc(table.n, sizeof(double), "q_q25");
    double *q_q50 = fc_alloc(table.n, sizeof(double), "q_q50");
    double *q_q75 = fc_alloc(table.n, sizeof(double), "q_q75");
    double *q_q95 = fc_alloc(table.n, sizeof(double), "q_q95");
    double *spd_mean = fc_alloc(table.n, sizeof(double), "spd_mean");
    double *spd_q05 = fc_alloc(table.n, sizeof(double), "spd_q05");
    double *spd_q25 = fc_alloc(table.n, sizeof(double), "spd_q25");
    double *spd_q50 = fc_alloc(table.n, sizeof(double), "spd_q50");
    double *spd_q75 = fc_alloc(table.n, sizeof(double), "spd_q75");
    double *spd_q95 = fc_alloc(table.n, sizeof(double), "spd_q95");

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
        len_q25[i] = row->length_quantiles.q25;
        len_q50[i] = row->length_quantiles.q50;
        len_q75[i] = row->length_quantiles.q75;
        len_q95[i] = row->length_quantiles.q95;
        q_mean[i] = row->qscore_mean;
        q_q05[i] = row->qscore_quantiles.q05;
        q_q25[i] = row->qscore_quantiles.q25;
        q_q50[i] = row->qscore_quantiles.q50;
        q_q75[i] = row->qscore_quantiles.q75;
        q_q95[i] = row->qscore_quantiles.q95;
        spd_mean[i] = row->speed_mean;
        spd_q05[i] = row->speed_quantiles.q05;
        spd_q25[i] = row->speed_quantiles.q25;
        spd_q50[i] = row->speed_quantiles.q50;
        spd_q75[i] = row->speed_quantiles.q75;
        spd_q95[i] = row->speed_quantiles.q95;
    }

    // Serialize as columnar table
    fc_json_object_t *data = fc_json_object_set_object(time, "__data__");
    fc_json_array_t *shape = fc_json_object_set_array(data, "shape");
    fc_json_array_append_number(shape, (double)table.n);

    fc_json_array_t *columns = fc_json_object_set_array(data, "columns");
    fc_json_array_append_value(columns, fc_json_column_f64("time_bin_start_s", time_start, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("time_bin_end_s", time_end, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("read_count", read_count, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("base_count", base_count, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("cumulative_read_count", cum_read, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("cumulative_base_count", cum_base, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_mean", len_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q05", len_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q25", len_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q50", len_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q75", len_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("length_q95", len_q95, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_mean", q_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q05", q_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q25", q_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q50", q_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q75", q_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("quality_q95", q_q95, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_mean", spd_mean, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q05", spd_q05, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q25", spd_q25, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q50", spd_q50, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q75", spd_q75, table.n));
    fc_json_array_append_value(columns, fc_json_column_f64("translocation_speed_q95", spd_q95, table.n));

    free(time_start);
    free(time_end);
    free(read_count);
    free(base_count);
    free(cum_read);
    free(cum_base);
    free(len_mean);
    free(len_q05);
    free(len_q25);
    free(len_q50);
    free(len_q75);
    free(len_q95);
    free(q_mean);
    free(q_q05);
    free(q_q25);
    free(q_q50);
    free(q_q75);
    free(q_q95);
    free(spd_mean);
    free(spd_q05);
    free(spd_q25);
    free(spd_q50);
    free(spd_q75);
    free(spd_q95);

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

    fc_json_object_t *channels = fc_json_object_set_object(root, "channels");
    fc_json_object_set_string(channels, "device", summary_channel_device_name(table.highest_channel));

    fc_channel_summary_skip_counts skipped = fc_channel_summary_skipped_counts(model);
    fc_json_object_t *skip = fc_json_object_set_object(channels, "skipped");
    fc_json_object_set_number(skip, "missing_channel", (double)skipped.missing_channel);
    fc_json_object_set_number(skip, "invalid_channel", (double)skipped.invalid_channel);
    fc_json_object_set_number(skip, "missing_mux", (double)skipped.missing_mux);
    fc_json_object_set_number(skip, "invalid_mux", (double)skipped.invalid_mux);

    // Pack channel data as schema-style table (list of column objects)
    fc_json_object_t *data = fc_json_object_set_object(channels, "__data__");
    fc_json_array_t *shape = fc_json_object_set_array(data, "shape");
    fc_json_array_append_number(shape, (double)table.n);

    size_t *channel_arr = fc_alloc(table.n, sizeof(size_t), "channel array");
    size_t *mux_arr = fc_alloc(table.n, sizeof(size_t), "mux array");
    size_t *read_count_arr = fc_alloc(table.n, sizeof(size_t), "read_count array");
    size_t *base_count_arr = fc_alloc(table.n, sizeof(size_t), "base_count array");

    for (size_t i = 0; i < table.n; i++) {
        const fc_channel_summary_row *row = &table.rows[i];
        channel_arr[i] = row->channel;
        mux_arr[i] = row->mux;
        read_count_arr[i] = row->read_count;
        base_count_arr[i] = row->base_count;
    }

    fc_json_array_t *columns = fc_json_object_set_array(data, "columns");
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("channel", channel_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("mux", mux_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("read_count", read_count_arr, table.n));
    fc_json_array_append_value(columns, fc_json_column_uint_adaptive("base_count", base_count_arr, table.n));

    free(channel_arr);
    free(mux_arr);
    free(read_count_arr);
    free(base_count_arr);

    fc_channel_summary_table_destroy(&table);
}

static void summary_json_flagstats(
    fc_json_object_t *root,
    const summary_json_sink *sink,
    const summary_json_barcode *state)
{
    if (!sink->cfg.bam_input || NULL == state->flagstats) {
        return;
    }

    fc_json_array_t *rows = fc_json_object_set_array(root, "flagstats");
    if (sink->cfg.region_mode) {
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *obj = fc_json_get_object(value);
        fc_json_object_set_string(obj, "ref", sink->cfg.region_label ? sink->cfg.region_label : "region");
        const size_t *counts = state->flagstats->counts[0];
        fc_json_object_set_number(obj, "total", (double)counts[0]);
        fc_json_object_set_number(obj, "primary", (double)counts[1]);
        fc_json_object_set_number(obj, "secondary", (double)counts[2]);
        fc_json_object_set_number(obj, "supplementary", (double)counts[3]);
        fc_json_object_set_number(obj, "unmapped", (double)counts[4]);
        fc_json_object_set_number(obj, "qcfail", (double)counts[5]);
        fc_json_object_set_number(obj, "duplicate", (double)counts[6]);
        fc_json_object_set_number(obj, "duplex", (double)counts[7]);
        fc_json_object_set_number(obj, "duplex_forming", (double)counts[8]);
        fc_json_array_append_value(rows, value);
    } else {
        for (size_t i = 0; i < state->flagstats->n_refs; i++) {
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
            const size_t *counts = state->flagstats->counts[i];
            fc_json_object_set_string(obj, "ref", label);
            fc_json_object_set_number(obj, "total", (double)counts[0]);
            fc_json_object_set_number(obj, "primary", (double)counts[1]);
            fc_json_object_set_number(obj, "secondary", (double)counts[2]);
            fc_json_object_set_number(obj, "supplementary", (double)counts[3]);
            fc_json_object_set_number(obj, "unmapped", (double)counts[4]);
            fc_json_object_set_number(obj, "qcfail", (double)counts[5]);
            fc_json_object_set_number(obj, "duplicate", (double)counts[6]);
            fc_json_object_set_number(obj, "duplex", (double)counts[7]);
            fc_json_object_set_number(obj, "duplex_forming", (double)counts[8]);
            fc_json_array_append_value(rows, value);
        }
    }
    if (NULL != state->flagstats->unmapped) {
        fc_json_value_t *value = fc_json_object_value();
        fc_json_object_t *obj = fc_json_get_object(value);
        const size_t *counts = state->flagstats->unmapped;
        fc_json_object_set_string(obj, "ref", "*");
        fc_json_object_set_number(obj, "total", (double)counts[0]);
        fc_json_object_set_number(obj, "primary", (double)counts[1]);
        fc_json_object_set_number(obj, "secondary", (double)counts[2]);
        fc_json_object_set_number(obj, "supplementary", (double)counts[3]);
        fc_json_object_set_number(obj, "unmapped", (double)counts[4]);
        fc_json_object_set_number(obj, "qcfail", (double)counts[5]);
        fc_json_object_set_number(obj, "duplicate", (double)counts[6]);
        fc_json_object_set_number(obj, "duplex", (double)counts[7]);
        fc_json_object_set_number(obj, "duplex_forming", (double)counts[8]);
        fc_json_array_append_value(rows, value);
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
    const read_stats *length_stats = summary_read_length_stats(engine, sink->cfg.unmapped_override);
    const read_stats *quality_stats = summary_quality_stats(engine, sink->cfg.bam_input, sink->cfg.unmapped_override);

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
    fc_json_object_set_number(summary, "n50_bp", (double)row.n50_bp);
    fc_json_object_set_number(summary, "l50", (double)row.l50);

    fc_json_object_t *sample = fc_json_object_set_object(root, "sample");
    summary_json_set_string_or_null(sample, "sample_id", row.sample_id);
    summary_json_set_string_or_null(sample, "input_type", row.input_type);
    summary_json_set_string_or_null(sample, "barcode_id", row.barcode_id);
    summary_json_set_string_or_null(sample, "barcode_alias", row.barcode_alias);

    summary_json_files(root, state);
    summary_json_counts(root, &row, state);
    summary_json_channels(root, sink->cfg.channel_summary, barcode);

    fc_json_object_t *distributions = fc_json_object_set_object(root, "distributions");
    summary_json_length_distribution(distributions, length_stats);
    summary_json_quality_distribution(distributions, quality_stats);

    summary_json_time(root, sink->cfg.time_bin, barcode);

    fc_json_object_t *heatmaps = fc_json_object_set_object(root, "heatmaps");
    summary_json_length_heatmap(heatmaps, sink->cfg.length_qscore, barcode);
    summary_json_speed_heatmap(heatmaps, sink->cfg.speed_qscore, barcode);
    summary_json_flagstats(root, sink, state);

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

    size_t barcode = (sink->cfg.use_barcodes && read_ctx) ? read_ctx->barcode : 0;
    if (barcode >= sink->cfg.barcode_count) {
        return;
    }
    summary_json_barcode *state = &sink->barcodes[barcode];
    summary_barcode_add_read(state, read);

    if (!sink->cfg.bam_input) {
        return;
    }
    summary_flag_stats *flagstats = summary_ensure_flagstats(sink, barcode);
    if (NULL == flagstats) {
        return;
    }

    bool is_unmapped = fc_is_unmapped(read->bam);
    if (!is_unmapped) {
        if (0 == flagstats->n_refs || NULL == flagstats->counts) {
            return;
        }
        size_t tid = (size_t)read->bam->core.tid;
        if (!sink->cfg.region_mode && tid >= flagstats->n_refs) {
            return;
        }
        size_t *counts = sink->cfg.region_mode ? flagstats->counts[0] : flagstats->counts[tid];
        summary_flagstats_add_counts(read->bam, counts, read->meta->duplex_code);
    }
    if (NULL != flagstats->unmapped && is_unmapped) {
        summary_flagstats_add_counts(read->bam, flagstats->unmapped, read->meta->duplex_code);
    }
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
    state->n_refs = state->cfg.region_mode ? 1 : (state->cfg.hdr ? (size_t)state->cfg.hdr->n_targets : 0);

    sink.ctx = state;
    sink.begin_file = summary_json_sink_begin_file;
    sink.add_read = summary_json_sink_add_read;
    sink.end_file = summary_json_sink_end_file;
    sink.finish = summary_json_sink_finish;
    sink.destroy = summary_json_sink_destroy;
    return sink;
}
