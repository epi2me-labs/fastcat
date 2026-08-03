#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "core/fc_time_bin.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_kll.h"

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
    int64_t first_bin;
    int64_t last_bin;
    fc_time_bin_accum *bins;
    size_t n_bins;
    size_t skipped_missing;
    size_t skipped_invalid;
    size_t skipped_duration_missing;
    size_t skipped_duration_invalid;
} fc_time_bin_barcode_state;

struct fc_time_bin {
    size_t barcode_count;
    double bin_size_s;
    size_t kll_k;
    bool origin_set;
    double origin_s;
    size_t skipped_missing;
    size_t skipped_invalid;
    size_t skipped_duration_missing;
    size_t skipped_duration_invalid;
    fc_time_bin_barcode_state **barcode_states;
};

bool fc_time_bin_parse_timestamp(const char *text, double *out)
{
    if (NULL == text || NULL == out || text[0] == '\0') {
        return false;
    }

    errno = 0;
    char *end = NULL;
    double value = strtod(text, &end);
    if (errno == 0 && end != text && *end == '\0' && isfinite(value)) {
        *out = value;
        return true;
    }

    struct tm tm = { 0 };
    tm.tm_isdst = 0;
    if (NULL == strptime(text, "%Y-%m-%dT%H:%M:%S", &tm)) {
        return false;
    }

    time_t epoch = timegm(&tm);
    if (epoch == (time_t)-1) {
        return false;
    }

    double epoch_s = (double)epoch;
    if (!isfinite(epoch_s)) {
        return false;
    }
    *out = epoch_s;
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

static void time_bin_destroy_barcode_state(fc_time_bin_barcode_state *state)
{
    if (NULL == state) {
        return;
    }
    for (size_t i = 0; i < state->n_bins; i++) {
        time_bin_destroy_accum(&state->bins[i]);
    }
    free(state->bins);
    free(state);
}

static fc_time_bin_accum *time_bin_ensure_accum(
    fc_time_bin *model,
    fc_time_bin_barcode_state *state,
    int64_t bin_idx)
{
    if (state->n_bins == 0) {
        state->bins = fc_alloc(1, sizeof(fc_time_bin_accum), "time_bin bins");
        state->n_bins = 1;
        state->first_bin = bin_idx;
        state->last_bin = bin_idx;
    } else if (bin_idx < state->first_bin) {
        size_t add = (size_t)(state->first_bin - bin_idx);
        size_t new_n = state->n_bins + add;
        fc_time_bin_accum *bins = fc_alloc(new_n, sizeof(*bins), "time_bin bins prepend");
        memcpy(bins + add, state->bins, state->n_bins * sizeof(*bins));
        free(state->bins);
        state->bins = bins;
        state->n_bins = new_n;
        state->first_bin = bin_idx;
    } else if (bin_idx > state->last_bin) {
        size_t old_n = state->n_bins;
        size_t add = (size_t)(bin_idx - state->last_bin);
        size_t new_n = old_n + add;
        state->bins = fc_recalloc(
            state->bins,
            old_n,
            new_n,
            sizeof(*state->bins),
            "time_bin bins append");
        state->n_bins = new_n;
        state->last_bin = bin_idx;
    }

    size_t offset = (size_t)(bin_idx - state->first_bin);
    fc_time_bin_accum *bin = &state->bins[offset];
    if (NULL == bin->length) {
        bin->length = fc_kll_create(model->kll_k);
    }
    if (NULL == bin->qscore) {
        bin->qscore = fc_kll_create(model->kll_k);
    }
    if (NULL == bin->speed) {
        bin->speed = fc_kll_create(model->kll_k);
    }
    return bin;
}

static fc_time_bin_accum *time_bin_find_accum(
    fc_time_bin *model,
    fc_time_bin_barcode_state *state,
    double timestamp_s)
{
    if (!model || !state || !model->origin_set || !isfinite(timestamp_s)) {
        return NULL;
    }
    double rel = (timestamp_s - model->origin_s) / model->bin_size_s;
    int64_t bin_idx = (int64_t)floor(rel);
    if (bin_idx < state->first_bin || bin_idx > state->last_bin) {
        return NULL;
    }
    return &state->bins[(size_t)(bin_idx - state->first_bin)];
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

    double timestamp_s = 0.0;
    if (!fc_time_bin_parse_timestamp(timestamp, &timestamp_s)) {
        state->skipped_invalid++;
        model->skipped_invalid++;
        return;
    }
    if (!isfinite(timestamp_s)) {
        return;
    }

    if (!model->origin_set) {
        model->origin_s = timestamp_s;
        model->origin_set = true;
    }

    double rel = (timestamp_s - model->origin_s) / model->bin_size_s;
    int64_t bin_idx = (int64_t)floor(rel);
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
    double timestamp_s = 0.0;
    if (!timestamp || timestamp[0] == '\0' || !fc_time_bin_parse_timestamp(timestamp, &timestamp_s)) {
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
    for (size_t i = 0; i < state->n_bins; ++i) {
        if (state->bins[i].accuracy_read_count > 0) {
            return true;
        }
    }
    return false;
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

static void time_bin_set_empty_distribution(fc_time_bin_row *row)
{
    row->length_mean = NAN;
    row->length_quantiles = time_bin_nan_quantiles();
    row->qscore_mean = NAN;
    row->qscore_quantiles = time_bin_nan_quantiles();
    row->speed_mean = NAN;
    row->speed_quantiles = time_bin_nan_quantiles();
    row->accuracy_mean = NAN;
    row->accuracy_quantiles = time_bin_nan_quantiles();
}

int fc_time_bin_materialize(const fc_time_bin *model, size_t barcode, fc_time_bin_table *out)
{
    if (NULL == model || NULL == out || barcode >= model->barcode_count) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    const fc_time_bin_barcode_state *state = model->barcode_states[barcode];
    if (NULL == state || !state->has_data || state->n_bins == 0) {
        return 0;
    }

    out->n = state->n_bins;
    out->rows = fc_alloc(out->n, sizeof(*out->rows), "time_bin rows");

    size_t cumulative_reads = 0;
    size_t cumulative_bases = 0;
    for (size_t i = 0; i < out->n; i++) {
        const fc_time_bin_accum *bin = &state->bins[i];
        fc_time_bin_row *row = &out->rows[i];

        row->time_bin_start_s = (double)i * model->bin_size_s;
        row->time_bin_end_s = row->time_bin_start_s + model->bin_size_s;
        row->read_count = bin->read_count;
        row->base_count = bin->base_count;
        cumulative_reads += bin->read_count;
        cumulative_bases += bin->base_count;
        row->cumulative_read_count = cumulative_reads;
        row->cumulative_base_count = cumulative_bases;

        if (bin->read_count == 0) {
            time_bin_set_empty_distribution(row);
            continue;
        }

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
    if (NULL == state || !state->has_data || state->n_bins == 0) {
        return 0;
    }
    if (max_rows == 0 || state->n_bins <= max_rows) {
        return fc_time_bin_materialize(model, barcode, out);
    }

    out->n = max_rows;
    out->rows = fc_alloc(out->n, sizeof(*out->rows), "time_bin merged rows");

    size_t cumulative_reads = 0;
    size_t cumulative_bases = 0;
    for (size_t i = 0; i < out->n; i++) {
        size_t start = i * state->n_bins / out->n;
        size_t end = (i + 1) * state->n_bins / out->n;
        fc_time_bin_row *row = &out->rows[i];

        row->time_bin_start_s = (double)start * model->bin_size_s;
        row->time_bin_end_s = (double)end * model->bin_size_s;

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
            const fc_time_bin_accum *bin = &state->bins[j];
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

        if (row->read_count == 0) {
            time_bin_set_empty_distribution(row);
        } else {
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
        }

        fc_kll_destroy(length);
        fc_kll_destroy(qscore);
        fc_kll_destroy(speed);
        fc_kll_destroy(accuracy);
    }

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
