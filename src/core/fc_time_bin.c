#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "core/fc_time_bin.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_kll.h"
#include "core/fc_time.h"
#include "core/kh_int64_map.h"

typedef struct {
    size_t read_count;
    size_t base_count;
    double sum_length;
    double sum_qscore;
    double sum_speed;
    size_t speed_count;
    fc_kll_sketch *length;
    fc_kll_sketch *qscore;
    fc_kll_sketch *speed;
    size_t accuracy_read_count;
    double sum_accuracy;
    fc_kll_sketch *accuracy;
} fc_time_bin_accum;

typedef struct {
    bool seen;
    bool has_data;
    kh_int64_map *bins;
    size_t skipped_missing;
    size_t skipped_invalid;
    size_t skipped_duration_missing;
    size_t skipped_duration_invalid;
} fc_time_bin_barcode_state;

struct fc_time_bin {
    size_t barcode_count;
    double bin_size_s;
    size_t kll_k;
    size_t skipped_missing;
    size_t skipped_invalid;
    size_t skipped_duration_missing;
    size_t skipped_duration_invalid;
    fc_time_bin_barcode_state **barcode_states;
};

typedef struct {
    int64_t index;
    const fc_time_bin_accum *accum;
} fc_time_bin_entry;

static bool time_bin_index(double timestamp_s, double bin_size_s, int64_t *out)
{
    // Absolute keys make arbitrary arrival order and wide timestamp gaps cheap.
    double index = floor(timestamp_s / bin_size_s);
    if (!isfinite(index) || index < (double)INT64_MIN || index >= -(double)INT64_MIN) {
        return false;
    }
    *out = (int64_t)index;
    return true;
}

static fc_time_bin_barcode_state *time_bin_ensure_barcode_state(fc_time_bin *model, size_t barcode)
{
    if (NULL == model || barcode >= model->barcode_count) {
        return NULL;
    }
    if (NULL == model->barcode_states[barcode]) {
        model->barcode_states[barcode] = fc_alloc(
            1,
            sizeof(fc_time_bin_barcode_state),
            "time_bin barcode state");
    }
    return model->barcode_states[barcode];
}

static void time_bin_destroy_accum(fc_time_bin_accum *bin)
{
    if (NULL == bin) {
        return;
    }
    fc_kll_destroy(bin->length);
    fc_kll_destroy(bin->qscore);
    fc_kll_destroy(bin->speed);
    fc_kll_destroy(bin->accuracy);
}

static void time_bin_destroy_accum_value(int64_t index, void *value, void *context)
{
    (void)index;
    (void)context;
    fc_time_bin_accum *bin = value;
    time_bin_destroy_accum(bin);
    free(bin);
}

static void time_bin_destroy_barcode_state(fc_time_bin_barcode_state *state)
{
    if (NULL == state) {
        return;
    }
    kh_int64_map_destroy(state->bins, time_bin_destroy_accum_value, NULL);
    free(state);
}

static fc_time_bin_accum *time_bin_ensure_accum(
    fc_time_bin *model,
    fc_time_bin_barcode_state *state,
    int64_t bin_idx)
{
    if (!state->bins) {
        state->bins = kh_int64_map_create();
        if (!state->bins) {
            return NULL;
        }
    }

    fc_time_bin_accum *existing = kh_int64_map_get(state->bins, bin_idx);
    if (existing) {
        return existing;
    }

    fc_time_bin_accum *bin = fc_alloc(1, sizeof(*bin), "time_bin accumulator");
    bin->length = fc_kll_create(model->kll_k);
    bin->qscore = fc_kll_create(model->kll_k);
    bin->speed = fc_kll_create(model->kll_k);
    if (!kh_int64_map_insert(state->bins, bin_idx, bin)) {
        time_bin_destroy_accum(bin);
        free(bin);
        return NULL;
    }
    return bin;
}

static fc_time_bin_accum *time_bin_find_accum(
    fc_time_bin *model,
    fc_time_bin_barcode_state *state,
    double timestamp_s)
{
    if (!model || !state || !state->bins || !isfinite(timestamp_s)) {
        return NULL;
    }
    int64_t bin_idx = 0;
    if (!time_bin_index(timestamp_s, model->bin_size_s, &bin_idx)) {
        return NULL;
    }
    return kh_int64_map_get(state->bins, bin_idx);
}

static int time_bin_entry_compare(const void *left, const void *right)
{
    const fc_time_bin_entry *a = left;
    const fc_time_bin_entry *b = right;
    return (a->index > b->index) - (a->index < b->index);
}

typedef struct {
    fc_time_bin_entry *entries;
    size_t out;
} fc_time_bin_entry_collector;

static void time_bin_collect_entry(int64_t index, void *value, void *context)
{
    fc_time_bin_entry_collector *collector = context;
    collector->entries[collector->out++] = (fc_time_bin_entry) {
        .index = index,
        .accum = value
    };
}

static void time_bin_note_accuracy(int64_t index, void *value, void *context)
{
    (void)index;
    const fc_time_bin_accum *bin = value;
    bool *seen = context;
    *seen = *seen || bin->accuracy_read_count > 0;
}

static fc_time_bin_entry *time_bin_sorted_entries(
    const fc_time_bin_barcode_state *state, size_t *n_entries)
{
    *n_entries = 0;
    if (!state || !state->bins || kh_int64_map_size(state->bins) == 0) {
        return NULL;
    }

    size_t n = kh_int64_map_size(state->bins);
    fc_time_bin_entry *entries = fc_alloc(n, sizeof(*entries), "time_bin entries");
    fc_time_bin_entry_collector collector = { .entries = entries };
    kh_int64_map_foreach(state->bins, time_bin_collect_entry, &collector);
    qsort(entries, collector.out, sizeof(*entries), time_bin_entry_compare);
    *n_entries = collector.out;
    return entries;
}

static double time_bin_elapsed_s(int64_t index, int64_t origin, double bin_size_s)
{
    return ((double)index - (double)origin) * bin_size_s;
}

fc_time_bin *fc_time_bin_create(size_t barcode_count, double bin_size_s, size_t kll_k)
{
    if (!isfinite(bin_size_s) || bin_size_s <= 0.0 || kll_k < 8) {
        return NULL;
    }

    fc_time_bin *model = fc_alloc(1, sizeof(*model), "time_bin model");
    model->barcode_count = barcode_count ? barcode_count : 1;
    model->bin_size_s = bin_size_s;
    model->kll_k = kll_k;
    model->barcode_states = fc_alloc(
        model->barcode_count,
        sizeof(fc_time_bin_barcode_state *),
        "time_bin barcode state array");
    return model;
}

void fc_time_bin_destroy(fc_time_bin *model)
{
    if (NULL == model) {
        return;
    }
    for (size_t i = 0; i < model->barcode_count; i++) {
        time_bin_destroy_barcode_state(model->barcode_states[i]);
    }
    free(model->barcode_states);
    free(model);
}

void fc_time_bin_add_read(fc_time_bin *model, size_t barcode, const fc_read *read)
{
    if (NULL == model || NULL == read || NULL == read->bam || NULL == read->meta
        || barcode >= model->barcode_count) {
        return;
    }
    if (!fc_is_canonical_read_record(read->bam)) {
        return;
    }

    fc_time_bin_barcode_state *state = time_bin_ensure_barcode_state(model, barcode);
    if (NULL == state) {
        return;
    }
    state->seen = true;

    const char *timestamp = read->meta->start_time;
    if (NULL == timestamp || timestamp[0] == '\0') {
        state->skipped_missing++;
        model->skipped_missing++;
        return;
    }

    double timestamp_s = read->meta->start_time_s;
    bool timestamp_valid = read->meta->start_time_valid;
    if (!read->meta->start_time_parsed) {
        timestamp_valid = fc_time_parse_timestamp(timestamp, &timestamp_s);
    }
    if (!timestamp_valid) {
        state->skipped_invalid++;
        model->skipped_invalid++;
        return;
    }
    if (!isfinite(timestamp_s)) {
        return;
    }

    int64_t bin_idx = 0;
    if (!time_bin_index(timestamp_s, model->bin_size_s, &bin_idx)) {
        state->skipped_invalid++;
        model->skipped_invalid++;
        return;
    }
    fc_time_bin_accum *bin = time_bin_ensure_accum(model, state, bin_idx);
    if (NULL == bin) {
        return;
    }

    double length = (double)read->bam->core.l_qseq;
    double mean_q = read->meta->mean_q;
    double duration_s = read->meta->duration_s;
    if (!isfinite(mean_q)) {
        mean_q = 0.0;
    }

    bin->read_count++;
    bin->base_count += (size_t)length;
    bin->sum_length += length;
    bin->sum_qscore += mean_q;
    fc_kll_add(bin->length, length);
    fc_kll_add(bin->qscore, mean_q);

    if (!isfinite(duration_s)) {
        state->skipped_duration_missing++;
        model->skipped_duration_missing++;
    } else if (duration_s <= 0.0) {
        state->skipped_duration_invalid++;
        model->skipped_duration_invalid++;
    } else {
        double speed = length / duration_s;
        if (!isfinite(speed)) {
            state->skipped_duration_invalid++;
            model->skipped_duration_invalid++;
        } else {
            bin->sum_speed += speed;
            bin->speed_count++;
            fc_kll_add(bin->speed, speed);
        }
    }
    state->has_data = true;
}

void fc_time_bin_add_accuracy(fc_time_bin *model, size_t barcode, const fc_read *read, float accuracy)
{
    if (!model || !read || !read->meta || !isfinite(accuracy) || barcode >= model->barcode_count) {
        return;
    }
    if (!read->bam || !fc_is_canonical_read_record(read->bam)) {
        return;
    }
    const char *timestamp = read->meta->start_time;
    if (!timestamp || timestamp[0] == '\0') {
        return;
    }
    double timestamp_s = read->meta->start_time_s;
    bool timestamp_valid = read->meta->start_time_valid;
    if (!read->meta->start_time_parsed) {
        timestamp_valid = fc_time_parse_timestamp(timestamp, &timestamp_s);
    }
    if (!timestamp_valid) {
        return;
    }
    int64_t bin_idx = 0;
    if (!time_bin_index(timestamp_s, model->bin_size_s, &bin_idx)) {
        return;
    }
    fc_time_bin_barcode_state *state = model->barcode_states[barcode];
    fc_time_bin_accum *bin = time_bin_find_accum(model, state, timestamp_s);
    if (!bin) {
        return;
    }
    if (!bin->accuracy) {
        bin->accuracy = fc_kll_create(model->kll_k);
    }
    double value = fmin(100.0, fmax(0.0, (double)accuracy));
    bin->accuracy_read_count++;
    bin->sum_accuracy += value;
    fc_kll_add(bin->accuracy, value);
}

bool fc_time_bin_seen(const fc_time_bin *model, size_t barcode)
{
    if (NULL == model || barcode >= model->barcode_count
        || NULL == model->barcode_states[barcode]) {
        return false;
    }
    return model->barcode_states[barcode]->seen;
}

bool fc_time_bin_accuracy_seen(const fc_time_bin *model, size_t barcode)
{
    if (NULL == model || barcode >= model->barcode_count) {
        return false;
    }
    const fc_time_bin_barcode_state *state = model->barcode_states[barcode];
    if (NULL == state) {
        return false;
    }
    if (!state->bins) {
        return false;
    }
    bool seen = false;
    kh_int64_map_foreach(state->bins, time_bin_note_accuracy, &seen);
    return seen;
}

fc_time_bin_skip_counts fc_time_bin_skipped_counts(const fc_time_bin *model)
{
    fc_time_bin_skip_counts counts = { 0 };
    if (NULL == model) {
        return counts;
    }
    counts.missing = model->skipped_missing;
    counts.invalid = model->skipped_invalid;
    counts.duration_missing = model->skipped_duration_missing;
    counts.duration_invalid = model->skipped_duration_invalid;
    return counts;
}

static fc_time_bin_quantiles time_bin_nan_quantiles(void)
{
    return (fc_time_bin_quantiles) {
        .q05 = NAN,
        .q10 = NAN,
        .q25 = NAN,
        .q50 = NAN,
        .q75 = NAN,
        .q90 = NAN,
        .q95 = NAN
    };
}

int fc_time_bin_materialize(const fc_time_bin *model, size_t barcode, fc_time_bin_table *out)
{
    if (NULL == model || NULL == out || barcode >= model->barcode_count) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    const fc_time_bin_barcode_state *state = model->barcode_states[barcode];
    if (NULL == state || !state->has_data) {
        return 0;
    }

    size_t n_entries = 0;
    fc_time_bin_entry *entries = time_bin_sorted_entries(state, &n_entries);
    if (!entries) {
        return 0;
    }

    out->n = n_entries;
    out->rows = fc_alloc(out->n, sizeof(*out->rows), "time_bin rows");

    size_t cumulative_reads = 0;
    size_t cumulative_bases = 0;
    for (size_t i = 0; i < out->n; i++) {
        const fc_time_bin_accum *bin = entries[i].accum;
        fc_time_bin_row *row = &out->rows[i];

        row->time_bin_start_s = time_bin_elapsed_s(
            entries[i].index, entries[0].index, model->bin_size_s);
        row->time_bin_end_s = row->time_bin_start_s + model->bin_size_s;
        row->read_count = bin->read_count;
        row->base_count = bin->base_count;
        cumulative_reads += bin->read_count;
        cumulative_bases += bin->base_count;
        row->cumulative_read_count = cumulative_reads;
        row->cumulative_base_count = cumulative_bases;

        row->length_mean = bin->sum_length / (double)bin->read_count;
        row->qscore_mean = bin->sum_qscore / (double)bin->read_count;
        row->length_quantiles = fc_kll_compute_quantile_septuple(bin->length);
        row->qscore_quantiles = fc_kll_compute_quantile_septuple(bin->qscore);
        if (bin->speed_count == 0) {
            row->speed_mean = NAN;
            row->speed_quantiles = time_bin_nan_quantiles();
        } else {
            row->speed_mean = bin->sum_speed / (double)bin->speed_count;
            row->speed_quantiles = fc_kll_compute_quantile_septuple(bin->speed);
        }
        row->accuracy_read_count = bin->accuracy_read_count;
        if (bin->accuracy_read_count > 0) {
            row->accuracy_mean = bin->sum_accuracy / (double)bin->accuracy_read_count;
            row->accuracy_quantiles = fc_kll_compute_quantile_septuple(bin->accuracy);
        } else {
            row->accuracy_mean = NAN;
            row->accuracy_quantiles = time_bin_nan_quantiles();
        }
    }

    free(entries);

    return 0;
}

int fc_time_bin_materialize_merged(
    const fc_time_bin *model,
    size_t barcode,
    size_t max_rows,
    fc_time_bin_table *out)
{
    if (NULL == model || NULL == out || barcode >= model->barcode_count) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    const fc_time_bin_barcode_state *state = model->barcode_states[barcode];
    if (NULL == state || !state->has_data) {
        return 0;
    }
    size_t n_entries = 0;
    fc_time_bin_entry *entries = time_bin_sorted_entries(state, &n_entries);
    if (!entries) {
        return 0;
    }
    if (max_rows == 0 || n_entries <= max_rows) {
        free(entries);
        return fc_time_bin_materialize(model, barcode, out);
    }

    out->n = max_rows;
    out->rows = fc_alloc(out->n, sizeof(*out->rows), "time_bin merged rows");

    size_t cumulative_reads = 0;
    size_t cumulative_bases = 0;
    for (size_t i = 0; i < out->n; i++) {
        size_t start = i * n_entries / out->n;
        size_t end = (i + 1) * n_entries / out->n;
        fc_time_bin_row *row = &out->rows[i];

        row->time_bin_start_s = time_bin_elapsed_s(
            entries[start].index, entries[0].index, model->bin_size_s);
        row->time_bin_end_s = time_bin_elapsed_s(
            entries[end - 1].index + 1, entries[0].index, model->bin_size_s);

        double sum_length = 0.0;
        double sum_qscore = 0.0;
        double sum_speed = 0.0;
        size_t speed_count = 0;
        double sum_accuracy = 0.0;
        size_t accuracy_read_count = 0;
        fc_kll_sketch *length = fc_kll_create(model->kll_k);
        fc_kll_sketch *qscore = fc_kll_create(model->kll_k);
        fc_kll_sketch *speed = fc_kll_create(model->kll_k);
        fc_kll_sketch *accuracy = fc_kll_create(model->kll_k);

        for (size_t j = start; j < end; j++) {
            const fc_time_bin_accum *bin = entries[j].accum;
            row->read_count += bin->read_count;
            row->base_count += bin->base_count;
            sum_length += bin->sum_length;
            sum_qscore += bin->sum_qscore;
            sum_speed += bin->sum_speed;
            speed_count += bin->speed_count;
            sum_accuracy += bin->sum_accuracy;
            accuracy_read_count += bin->accuracy_read_count;
            fc_kll_merge(length, bin->length);
            fc_kll_merge(qscore, bin->qscore);
            fc_kll_merge(speed, bin->speed);
            fc_kll_merge(accuracy, bin->accuracy);
        }

        cumulative_reads += row->read_count;
        cumulative_bases += row->base_count;
        row->cumulative_read_count = cumulative_reads;
        row->cumulative_base_count = cumulative_bases;

        row->length_mean = sum_length / (double)row->read_count;
        row->qscore_mean = sum_qscore / (double)row->read_count;
        row->length_quantiles = fc_kll_compute_quantile_septuple(length);
        row->qscore_quantiles = fc_kll_compute_quantile_septuple(qscore);
        if (speed_count == 0) {
            row->speed_mean = NAN;
            row->speed_quantiles = time_bin_nan_quantiles();
        } else {
            row->speed_mean = sum_speed / (double)speed_count;
            row->speed_quantiles = fc_kll_compute_quantile_septuple(speed);
        }
        row->accuracy_read_count = accuracy_read_count;
        if (accuracy_read_count == 0) {
            row->accuracy_mean = NAN;
            row->accuracy_quantiles = time_bin_nan_quantiles();
        } else {
            row->accuracy_mean = sum_accuracy / (double)accuracy_read_count;
            row->accuracy_quantiles = fc_kll_compute_quantile_septuple(accuracy);
        }

        fc_kll_destroy(length);
        fc_kll_destroy(qscore);
        fc_kll_destroy(speed);
        fc_kll_destroy(accuracy);
    }

    free(entries);

    return 0;
}

void fc_time_bin_table_destroy(fc_time_bin_table *table)
{
    if (NULL == table) {
        return;
    }
    free(table->rows);
    table->rows = NULL;
    table->n = 0;
}
