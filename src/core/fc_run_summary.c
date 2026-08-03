#include "core/fc_run_summary.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_stats_engine.h"
#include "core/fc_time_bin.h"
#include "core/kh_counter.h"

typedef struct {
    bool seen;
    kh_counter_t *sample_id;
    kh_counter_t *run_id;
    kh_counter_t *basecaller_model;
    kh_counter_t *modbase_models;
    kh_counter_t *flow_cell_id;
    kh_counter_t *barcode_alias;
    char *run_date;
    double min_time_s;
    double max_time_s;
    bool has_time;
    size_t read_count;
    size_t base_count;
    uint64_t filter_counts[NUM_FAILURE_CODES];
} fc_run_summary_state;

struct fc_run_summary {
    size_t barcode_count;
    fc_run_summary_state *states;
};

static char *summary_join_counter(const kh_counter_t *counter, const char *label)
{
    if (NULL == counter || 0 == kh_size(counter)) {
        return NULL;
    }

    size_t total_len = 0;
    size_t count = 0;
    for (khiter_t k = 0; k < kh_end(counter); ++k) {
        if (!kh_exist(counter, k)) {
            continue;
        }
        total_len += strlen(kh_key(counter, k));
        if (count > 0) {
            total_len += 1;
        }
        count += 1;
    }

    char *joined = fc_alloc(total_len + 1, sizeof(char), label);
    char *dst = joined;
    count = 0;
    for (khiter_t k = 0; k < kh_end(counter); ++k) {
        if (!kh_exist(counter, k)) {
            continue;
        }
        const char *value = kh_key(counter, k);
        size_t len = strlen(value);
        if (count > 0) {
            *dst++ = ';';
        }
        memcpy(dst, value, len);
        dst += len;
        count += 1;
    }
    *dst = '\0';
    return joined;
}

static char *summary_format_timestamp(const char *value)
{
    if (!value || value[0] == '\0') {
        return fc_strdup("", "run_summary timestamp");
    }
    if (strlen(value) >= 19
        && value[4] == '-'
        && value[7] == '-'
        && (value[10] == 'T' || value[10] == ' ')
        && value[13] == ':'
        && value[16] == ':') {
        return fc_asprintf(
            "run_summary timestamp",
            "%.10s %.8s",
            value,
            value + 11);
    }
    return fc_strdup(value, "run_summary timestamp");
}

static size_t length_bin_start_bp(const read_stats *stats, size_t index)
{
    if (!stats) {
        return 0;
    }
    if (stats->width == 0.0f && stats->edges) {
        return stats->edges[index];
    }
    return (size_t)llround((double)index * (double)stats->width);
}

static void compute_nxx_lxx(const read_stats *length_stats, fc_run_summary_row *out)
{
    const double thresholds[] = { 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95 };
    size_t *nxx[] = { &out->n5_bp, &out->n10_bp, &out->n25_bp, &out->n50_bp,
        &out->n75_bp, &out->n90_bp, &out->n95_bp };
    size_t *lxx[] = { &out->l5, &out->l10, &out->l25, &out->l50,
        &out->l75, &out->l90, &out->l95 };
    const size_t n_thresholds = sizeof(thresholds) / sizeof(thresholds[0]);
    for (size_t i = 0; i < n_thresholds; ++i) {
        *nxx[i] = 0;
        *lxx[i] = 0;
    }
    if (!length_stats || !length_stats->counts || length_stats->n == 0) {
        return;
    }

    double total_bases = 0.0;
    for (size_t i = 0; i < length_stats->n; ++i) {
        size_t bin_start = length_bin_start_bp(length_stats, i);
        total_bases += (double)bin_start * (double)length_stats->counts[i];
    }
    if (total_bases <= 0.0) {
        return;
    }

    double cumulative_bases = 0.0;
    size_t cumulative_reads = 0;
    size_t next_threshold = 0;
    for (size_t i = length_stats->n; i-- > 0 && next_threshold < n_thresholds;) {
        size_t count = length_stats->counts[i];
        if (count == 0) {
            continue;
        }
        size_t bin_start = length_bin_start_bp(length_stats, i);
        cumulative_bases += (double)bin_start * (double)count;
        cumulative_reads += count;
        while (next_threshold < n_thresholds && cumulative_bases >= total_bases * thresholds[next_threshold]) {
            *nxx[next_threshold] = bin_start;
            *lxx[next_threshold] = cumulative_reads;
            ++next_threshold;
        }
    }
}

static double summary_stats_total_count(const read_stats *stats)
{
    if (NULL == stats) {
        return 0.0;
    }
    double total = 0.0;
    for (size_t i = 0; i < stats->n; ++i) {
        total += (double)stats->counts[i];
    }
    return total;
}

static const char *summary_input_type(const fc_stats_engine *engine, bool bam_input)
{
    if (!bam_input) {
        return "fastq";
    }
    // if acc_stats is present and has data, we saw at least one aligned read
    if (NULL != engine && summary_stats_total_count(engine->acc_stats) > 0.0) {
        return "aligned bam";
    }
    return "unaligned bam";
}

fc_run_summary *fc_run_summary_create(size_t barcode_count)
{
    fc_run_summary *model = fc_alloc(1, sizeof(*model), "run_summary model");
    model->barcode_count = barcode_count ? barcode_count : 1;
    model->states = fc_alloc(model->barcode_count, sizeof(*model->states), "run_summary states");
    return model;
}

void fc_run_summary_destroy(fc_run_summary *model)
{
    if (!model) {
        return;
    }
    for (size_t i = 0; i < model->barcode_count; ++i) {
        kh_counter_destroy(model->states[i].run_id);
        kh_counter_destroy(model->states[i].sample_id);
        kh_counter_destroy(model->states[i].basecaller_model);
        kh_counter_destroy(model->states[i].modbase_models);
        kh_counter_destroy(model->states[i].flow_cell_id);
        kh_counter_destroy(model->states[i].barcode_alias);
        free(model->states[i].run_date);
    }
    free(model->states);
    free(model);
}

void fc_run_summary_row_destroy(fc_run_summary_row *row)
{
    if (NULL == row) {
        return;
    }
    free(row->sample_id);
    free(row->input_type);
    free(row->barcode_id);
    free(row->barcode_alias);
    free(row->run_id);
    free(row->basecaller_model);
    free(row->modbase_models);
    free(row->flow_cell_id);
    free(row->run_date);
    memset(row, 0, sizeof(*row));
}

void fc_run_summary_add_read(fc_run_summary *model, size_t barcode, const fc_read *read)
{
    if (!model || !read || !read->bam || barcode >= model->barcode_count) {
        return;
    }
    if (!fc_is_canonical_read_record(read->bam)) {
        return;
    }

    fc_run_summary_state *state = &model->states[barcode];
    state->seen = true;
    state->read_count += 1;
    state->base_count += read->bam->core.l_qseq;

    if (NULL != read->meta->sample_id && read->meta->sample_id[0] != '\0') {
        if (NULL == state->sample_id) {
            state->sample_id = kh_counter_init();
        }
        kh_counter_increment(state->sample_id, read->meta->sample_id);
    }
    if (NULL != read->meta->runid && read->meta->runid[0] != '\0') {
        if (NULL == state->run_id) {
            state->run_id = kh_counter_init();
        }
        kh_counter_increment(state->run_id, read->meta->runid);
    }
    if (NULL != read->meta->basecaller && read->meta->basecaller[0] != '\0') {
        if (NULL == state->basecaller_model) {
            state->basecaller_model = kh_counter_init();
        }
        kh_counter_increment(state->basecaller_model, read->meta->basecaller);
    }
    if (NULL != read->meta->modcallers && read->meta->modcallers[0] != '\0') {
        if (NULL == state->modbase_models) {
            state->modbase_models = kh_counter_init();
        }
        kh_counter_increment(state->modbase_models, read->meta->modcallers);
    }
    if (NULL != read->meta->flow_cell_id && read->meta->flow_cell_id[0] != '\0') {
        if (NULL == state->flow_cell_id) {
            state->flow_cell_id = kh_counter_init();
        }
        kh_counter_increment(state->flow_cell_id, read->meta->flow_cell_id);
    }
    if (NULL != read->meta->barcode_alias && read->meta->barcode_alias[0] != '\0') {
        if (NULL == state->barcode_alias) {
            state->barcode_alias = kh_counter_init();
        }
        kh_counter_increment(state->barcode_alias, read->meta->barcode_alias);
    }

    if (NULL != read->meta->start_time && read->meta->start_time[0] != '\0') {
        double time_s = 0.0;
        if (fc_time_bin_parse_timestamp(read->meta->start_time, &time_s)) {
            if (!state->has_time || time_s < state->min_time_s) {
                state->min_time_s = time_s;
                free(state->run_date);
                state->run_date = summary_format_timestamp(read->meta->start_time);
            }
            if (!state->has_time || time_s > state->max_time_s) {
                state->max_time_s = time_s;
            }
            state->has_time = true;
        }
    }
}

void fc_run_summary_add_filter_event(fc_run_summary *model, size_t barcode, failure_code code)
{
    if (!model || barcode >= model->barcode_count || code >= NUM_FAILURE_CODES) {
        return;
    }
    fc_run_summary_state *state = &model->states[barcode];
    state->seen = true;
    state->filter_counts[code]++;
}

bool fc_run_summary_seen(const fc_run_summary *model, size_t barcode)
{
    return NULL != model && barcode < model->barcode_count && model->states[barcode].seen;
}

int fc_run_summary_materialize(
    const fc_run_summary *model,
    size_t barcode,
    const fc_stats_engine *engine,
    bool bam_input,
    const read_stats *length_stats,
    const char *sample_id,
    const char *barcode_id,
    fc_run_summary_row *out)
{
    if (!model || !out || barcode >= model->barcode_count) {
        return -1;
    }

    fc_run_summary_state *state = &model->states[barcode];
    memset(out, 0, sizeof(*out));
    out->sample_id = NULL != state->sample_id && kh_size(state->sample_id) > 0
        ? summary_join_counter(state->sample_id, "run_summary sample id")
        : fc_strdup(sample_id ? sample_id : "SAMPLE", "run_summary sample fallback");
    out->input_type = fc_strdup(summary_input_type(engine, bam_input), "run_summary input type");
    out->barcode_id = fc_strdup(barcode_id ? barcode_id : "", "run_summary barcode id");
    out->barcode_alias = summary_join_counter(state->barcode_alias, "run_summary barcode alias");
    if (NULL == out->barcode_alias) {
        out->barcode_alias = fc_strdup("", "run_summary barcode alias empty");
    }
    out->run_id = summary_join_counter(state->run_id, "run_summary run id");
    out->basecaller_model = summary_join_counter(state->basecaller_model, "run_summary basecaller");
    out->modbase_models = summary_join_counter(state->modbase_models, "run_summary modbase models");
    out->flow_cell_id = summary_join_counter(state->flow_cell_id, "run_summary flow cell id");
    out->run_date = NULL != state->flow_cell_id && kh_size(state->flow_cell_id) > 1
        ? fc_strdup("(not meaningful)", "run_summary run date")
        : fc_strdup(state->run_date ? state->run_date : "", "run_summary run date");
    if (NULL == out->run_id) {
        out->run_id = fc_strdup("", "run_summary run id empty");
    }
    if (NULL == out->basecaller_model) {
        out->basecaller_model = fc_strdup("", "run_summary basecaller empty");
    }
    if (NULL == out->modbase_models) {
        out->modbase_models = fc_strdup("", "run_summary modbase empty");
    }
    if (NULL == out->flow_cell_id) {
        out->flow_cell_id = fc_strdup("", "run_summary flow cell id empty");
    }
    out->run_duration_s = NULL != state->flow_cell_id && kh_size(state->flow_cell_id) > 1
        ? -1.0
        : (state->has_time ? fmax(0.0, state->max_time_s - state->min_time_s) : NAN);
    out->read_count = state->read_count;
    out->base_count = state->base_count;
    out->records_seen = state->filter_counts[R_RECORD_SEEN];
    out->records_passed_filters = state->filter_counts[R_RECORD_OK];
    out->records_too_long = state->filter_counts[R_TOO_LONG];
    out->records_too_short = state->filter_counts[R_TOO_SHORT];
    out->records_low_quality = state->filter_counts[R_LOW_QUALITY];
    out->records_dust_masked = state->filter_counts[R_DUST_MASKED];
    compute_nxx_lxx(length_stats, out);
    return 0;
}
