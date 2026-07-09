#ifndef FASTCAT_FC_RUN_SUMMARY_H
#define FASTCAT_FC_RUN_SUMMARY_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_read.h"
#include "core/fc_stats.h"
#include "core/fc_stats_engine.h"

typedef struct fc_run_summary fc_run_summary;

typedef struct {
    const char *sample_id;
    const char *input_type;
    const char *barcode_id;
    const char *barcode_alias;
    const char *run_id;
    const char *basecaller_model;
    const char *flow_cell_id;
    const char *run_date;
    double run_duration_s;
    size_t read_count;
    size_t base_count;
    size_t n50_bp;
    size_t l50;
} fc_run_summary_row;

fc_run_summary *fc_run_summary_create(size_t barcode_count);
void fc_run_summary_destroy(fc_run_summary *model);

void fc_run_summary_add_read(fc_run_summary *model, size_t barcode, const fc_read *read);

bool fc_run_summary_seen(const fc_run_summary *model, size_t barcode);

int fc_run_summary_materialize(
    const fc_run_summary *model,
    size_t barcode,
    const fc_stats_engine *engine,
    bool bam_input,
    const read_stats *length_stats,
    const char *sample_id,
    const char *barcode_id,
    fc_run_summary_row *out);

#endif
