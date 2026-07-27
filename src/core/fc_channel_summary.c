#include "core/fc_channel_summary.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"

typedef struct {
    bool seen;
    size_t highest_channel;
    size_t *read_counts;
    size_t *base_counts;
    size_t *quality_read_counts;
    double *quality_sums;
    size_t *concordance_read_counts;
    double *concordance_sums;
    fc_channel_summary_skip_counts skipped;
} fc_channel_summary_barcode_state;

struct fc_channel_summary {
    size_t barcode_count;
    fc_channel_summary_barcode_state **barcode_states;
    fc_channel_summary_skip_counts skipped;
};

static size_t channel_summary_cell_index(size_t channel, size_t mux)
{
    return ((channel - 1) * FC_CHANNEL_SUMMARY_MUX_SLOTS) + mux;
}

static fc_channel_summary_barcode_state *channel_summary_ensure_barcode_state(
    fc_channel_summary *model,
    size_t barcode)
{
    if (NULL == model || barcode >= model->barcode_count) {
        return NULL;
    }
    if (NULL == model->barcode_states[barcode]) {
        fc_channel_summary_barcode_state *state = fc_alloc(
            1, sizeof(*state), "channel_summary barcode state");
        size_t n_cells = FC_CHANNEL_SUMMARY_MAX_CHANNEL * FC_CHANNEL_SUMMARY_MUX_SLOTS;
        state->read_counts = fc_alloc(n_cells, sizeof(size_t), "channel_summary read counts");
        state->base_counts = fc_alloc(n_cells, sizeof(size_t), "channel_summary base counts");
        state->quality_read_counts = fc_alloc(n_cells, sizeof(size_t), "channel_summary quality counts");
        state->quality_sums = fc_alloc(n_cells, sizeof(double), "channel_summary quality sums");
        state->concordance_read_counts = fc_alloc(n_cells, sizeof(size_t), "channel_summary concordance counts");
        state->concordance_sums = fc_alloc(n_cells, sizeof(double), "channel_summary concordance sums");
        model->barcode_states[barcode] = state;
    }
    return model->barcode_states[barcode];
}

static void channel_summary_destroy_barcode_state(fc_channel_summary_barcode_state *state)
{
    if (NULL == state) {
        return;
    }
    free(state->read_counts);
    free(state->base_counts);
    free(state->quality_read_counts);
    free(state->quality_sums);
    free(state->concordance_read_counts);
    free(state->concordance_sums);
    free(state);
}

static bool channel_summary_concordance_qscore(const fc_read *read, double *out)
{
    if (!fc_is_primary_qc_pass(read->bam) || read->meta->edit_distance < 0) {
        return false;
    }

    size_t match = 0;
    size_t ins = 0;
    size_t del = 0;
    uint32_t *cigar = bam_get_cigar(read->bam);
    for (uint32_t i = 0; i < read->bam->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        size_t length = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            match += length;
        } else if (op == BAM_CINS) {
            ins += length;
        } else if (op == BAM_CDEL) {
            del += length;
        }
    }

    size_t aligned_length = match + ins + del;
    size_t edit_distance = (size_t)read->meta->edit_distance;
    if (aligned_length == 0 || edit_distance > aligned_length) {
        return false;
    }
    double error_rate = (double)edit_distance / (double)aligned_length;
    *out = error_rate > 0.0 ? -10.0 * log10(error_rate) : 60.0;
    *out = fmin(60.0, *out);
    return isfinite(*out);
}

fc_channel_summary *fc_channel_summary_create(size_t barcode_count)
{
    fc_channel_summary *model = fc_alloc(1, sizeof(*model), "channel_summary");
    model->barcode_count = barcode_count ? barcode_count : 1;
    model->barcode_states = fc_alloc(
        model->barcode_count,
        sizeof(fc_channel_summary_barcode_state *),
        "channel_summary barcode states");
    return model;
}

void fc_channel_summary_destroy(fc_channel_summary *model)
{
    if (NULL == model) {
        return;
    }
    for (size_t i = 0; i < model->barcode_count; i++) {
        channel_summary_destroy_barcode_state(model->barcode_states[i]);
    }
    free(model->barcode_states);
    free(model);
}

void fc_channel_summary_add_read(fc_channel_summary *model, size_t barcode, const fc_read *read)
{
    if (NULL == model || NULL == read || NULL == read->bam || NULL == read->meta
        || barcode >= model->barcode_count) {
        return;
    }

    fc_channel_summary_barcode_state *state = channel_summary_ensure_barcode_state(model, barcode);
    if (NULL == state) {
        return;
    }
    state->seen = true;

    size_t channel = read->meta->channel;
    size_t mux = read->meta->mux;

    if (channel == 0) {
        state->skipped.missing_channel++;
        model->skipped.missing_channel++;
        return;
    }
    if (channel > FC_CHANNEL_SUMMARY_MAX_CHANNEL) {
        state->skipped.invalid_channel++;
        model->skipped.invalid_channel++;
        return;
    }
    if (mux > FC_CHANNEL_SUMMARY_MAX_MUX) {
        state->skipped.invalid_mux++;
        model->skipped.invalid_mux++;
        return;
    }

    size_t idx = channel_summary_cell_index(channel, mux);
    state->read_counts[idx]++;
    state->base_counts[idx] += (size_t)read->bam->core.l_qseq;
    if (isfinite(read->meta->mean_q) && read->meta->mean_q >= 0.0f) {
        state->quality_read_counts[idx]++;
        state->quality_sums[idx] += read->meta->mean_q;
    }
    double concordance_qscore = 0.0;
    if (channel_summary_concordance_qscore(read, &concordance_qscore)) {
        state->concordance_read_counts[idx]++;
        state->concordance_sums[idx] += concordance_qscore;
    }
    if (channel > state->highest_channel) {
        state->highest_channel = channel;
    }
}

bool fc_channel_summary_seen(const fc_channel_summary *model, size_t barcode)
{
    if (NULL == model || barcode >= model->barcode_count
        || NULL == model->barcode_states[barcode]) {
        return false;
    }
    return model->barcode_states[barcode]->seen;
}

fc_channel_summary_skip_counts fc_channel_summary_skipped_counts(const fc_channel_summary *model)
{
    fc_channel_summary_skip_counts counts = { 0 };
    if (NULL == model) {
        return counts;
    }
    return model->skipped;
}

int fc_channel_summary_materialize_sparse(
    const fc_channel_summary *model,
    size_t barcode,
    fc_channel_summary_table *out)
{
    if (NULL == model || NULL == out || barcode >= model->barcode_count) {
        return -1;
    }

    memset(out, 0, sizeof(*out));
    const fc_channel_summary_barcode_state *state = model->barcode_states[barcode];
    if (NULL == state || state->highest_channel == 0) {
        return 0;
    }

    size_t n = 0;
    for (size_t channel = 1; channel <= state->highest_channel; channel++) {
        for (size_t mux = FC_CHANNEL_SUMMARY_UNKNOWN_MUX; mux <= FC_CHANNEL_SUMMARY_MAX_MUX; mux++) {
            size_t idx = channel_summary_cell_index(channel, mux);
            if (state->read_counts[idx] != 0 || state->base_counts[idx] != 0) {
                n++;
            }
        }
    }

    out->highest_channel = state->highest_channel;
    out->n = n;
    if (n == 0) {
        return 0;
    }

    out->rows = fc_alloc(n, sizeof(*out->rows), "channel_summary rows");
    size_t row = 0;
    for (size_t channel = 1; channel <= state->highest_channel; channel++) {
        for (size_t mux = FC_CHANNEL_SUMMARY_UNKNOWN_MUX; mux <= FC_CHANNEL_SUMMARY_MAX_MUX; mux++) {
            size_t idx = channel_summary_cell_index(channel, mux);
            size_t read_count = state->read_counts[idx];
            size_t base_count = state->base_counts[idx];
            if (read_count == 0 && base_count == 0) {
                continue;
            }
            out->rows[row++] = (fc_channel_summary_row) {
                .channel = channel,
                .mux = mux,
                .read_count = read_count,
                .base_count = base_count,
                .quality_read_count = state->quality_read_counts[idx],
                .mean_quality = state->quality_read_counts[idx] > 0
                    ? state->quality_sums[idx] / (double)state->quality_read_counts[idx]
                    : NAN,
                .concordance_read_count = state->concordance_read_counts[idx],
                .mean_concordance_qscore = state->concordance_read_counts[idx] > 0
                    ? state->concordance_sums[idx] / (double)state->concordance_read_counts[idx]
                    : NAN
            };
        }
    }
    return 0;
}

void fc_channel_summary_table_destroy(fc_channel_summary_table *table)
{
    if (NULL == table) {
        return;
    }
    free(table->rows);
    memset(table, 0, sizeof(*table));
}
