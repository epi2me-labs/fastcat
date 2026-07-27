#ifndef FASTCAT_FC_CHANNEL_SUMMARY_H
#define FASTCAT_FC_CHANNEL_SUMMARY_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_read.h"

#define FC_CHANNEL_SUMMARY_MAX_CHANNEL 3000
#define FC_CHANNEL_SUMMARY_UNKNOWN_MUX 0
#define FC_CHANNEL_SUMMARY_MAX_MUX 4
#define FC_CHANNEL_SUMMARY_MUX_SLOTS (FC_CHANNEL_SUMMARY_MAX_MUX + 1)

typedef struct fc_channel_summary fc_channel_summary;

typedef struct {
    size_t channel;
    size_t mux;
    size_t read_count;
    size_t base_count;
    size_t quality_read_count;
    double mean_quality;
    size_t concordance_read_count;
    double mean_concordance_qscore;
} fc_channel_summary_row;

typedef struct {
    fc_channel_summary_row *rows;
    size_t n;
    size_t highest_channel;
} fc_channel_summary_table;

typedef struct {
    size_t missing_channel;
    size_t invalid_channel;
    size_t missing_mux;
    size_t invalid_mux;
} fc_channel_summary_skip_counts;

fc_channel_summary *fc_channel_summary_create(size_t barcode_count);
void fc_channel_summary_destroy(fc_channel_summary *model);

void fc_channel_summary_add_read(fc_channel_summary *model, size_t barcode, const fc_read *read);

bool fc_channel_summary_seen(const fc_channel_summary *model, size_t barcode);
fc_channel_summary_skip_counts fc_channel_summary_skipped_counts(const fc_channel_summary *model);

int fc_channel_summary_materialize_sparse(
    const fc_channel_summary *model,
    size_t barcode,
    fc_channel_summary_table *out);
void fc_channel_summary_table_destroy(fc_channel_summary_table *table);

#endif
