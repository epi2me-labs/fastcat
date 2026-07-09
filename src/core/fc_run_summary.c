#include "core/fc_run_summary.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_stats_engine.h"
#include "core/fc_time_bin.h"

typedef struct {
    bool seen;
    char *sample_id;
    char *run_id;
    char *basecaller_model;
    char *flow_cell_id;
    char *barcode_alias;
    bool sample_id_multiple;
    bool run_id_multiple;
    bool basecaller_multiple;
    bool flow_cell_multiple;
    char *run_date;
    double min_time_s;
    double max_time_s;
    bool has_time;
    size_t read_count;
    size_t base_count;
} fc_run_summary_state;

struct fc_run_summary {
    size_t barcode_count;
    fc_run_summary_state *states;
};

static void summary_update_string(char **dst, bool *multiple, const char *src, const char *label)
{
    if (!src || src[0] == '\0') {
        return;
    }
    if (!*dst || (*dst)[0] == '\0') {
        *dst = fc_strdup(src, label);
        return;
    }
    if (strcmp(*dst, src) != 0) {
        *multiple = true;
    }
}

static void summary_update_first_string(char **dst, const char *src, const char *label)
{
    if (!src || src[0] == '\0' || (*dst && (*dst)[0] != '\0')) {
        return;
    }
    *dst = fc_strdup(src, label);
}

static const char *summary_display_string(const char *value, bool multiple)
{
    if (multiple) {
        return "(multiple)";
    }
    return value ? value : "";
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

static void compute_n50_l50(const read_stats *length_stats, size_t *n50_bp, size_t *l50)
{
    *n50_bp = 0;
    *l50 = 0;
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

    double target = total_bases * 0.5;
    double cumulative_bases = 0.0;
    size_t cumulative_reads = 0;
    for (size_t i = 0; i < length_stats->n; ++i) {
        size_t count = length_stats->counts[i];
        if (count == 0) {
            continue;
        }
        size_t bin_start = length_bin_start_bp(length_stats, i);
        cumulative_bases += (double)bin_start * (double)count;
        cumulative_reads += count;
        if (cumulative_bases >= target) {
            *n50_bp = bin_start;
            *l50 = cumulative_reads;
            return;
        }
    }
}

static double summary_stats_total_count(const read_stats *stats)
{
    if (!stats) {
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
    if (engine && summary_stats_total_count(engine->length_stats) > 0.0) {
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
        free(model->states[i].run_id);
        free(model->states[i].sample_id);
        free(model->states[i].basecaller_model);
        free(model->states[i].flow_cell_id);
        free(model->states[i].barcode_alias);
        free(model->states[i].run_date);
    }
    free(model->states);
    free(model);
}

void fc_run_summary_add_read(fc_run_summary *model, size_t barcode, const fc_read *read)
{
    if (!model || !read || !read->bam || barcode >= model->barcode_count) {
        return;
    }

    fc_run_summary_state *state = &model->states[barcode];
    state->seen = true;
    state->read_count += 1;
    state->base_count += read->bam->core.l_qseq;

    summary_update_string(&state->sample_id, &state->sample_id_multiple, read->meta->sample_id, "run_summary sample id");
    summary_update_string(&state->run_id, &state->run_id_multiple, read->meta->runid, "run_summary run id");
    summary_update_string(&state->basecaller_model, &state->basecaller_multiple, read->meta->basecaller, "run_summary basecaller");
    summary_update_string(&state->flow_cell_id, &state->flow_cell_multiple, read->meta->flow_cell_id, "run_summary flow cell id");
    summary_update_first_string(&state->barcode_alias, read->meta->barcode_alias, "run_summary barcode alias");

    if (read->meta->start_time && read->meta->start_time[0] != '\0') {
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

bool fc_run_summary_seen(const fc_run_summary *model, size_t barcode)
{
    return model && barcode < model->barcode_count && model->states[barcode].seen;
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

    const fc_run_summary_state *state = &model->states[barcode];
    memset(out, 0, sizeof(*out));
    out->sample_id = state->sample_id_multiple
        ? "(multiple)"
        : ((state->sample_id && state->sample_id[0] != '\0') ? state->sample_id : (sample_id ? sample_id : ""));
    out->input_type = summary_input_type(engine, bam_input);
    out->barcode_id = barcode_id ? barcode_id : "";
    out->barcode_alias = state->barcode_alias ? state->barcode_alias : "";
    out->run_id = summary_display_string(state->run_id, state->run_id_multiple);
    out->basecaller_model = summary_display_string(state->basecaller_model, state->basecaller_multiple);
    out->flow_cell_id = summary_display_string(state->flow_cell_id, state->flow_cell_multiple);
    out->run_date = state->flow_cell_multiple
        ? "(not meaningful)"
        : (state->run_date ? state->run_date : "");
    out->run_duration_s = state->flow_cell_multiple
        ? -1.0
        : (state->has_time ? fmax(0.0, state->max_time_s - state->min_time_s) : NAN);
    out->read_count = state->read_count;
    out->base_count = state->base_count;
    compute_n50_l50(length_stats, &out->n50_bp, &out->l50);
    return 0;
}
