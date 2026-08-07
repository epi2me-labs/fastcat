#ifndef FASTCAT_FC_RUN_SUMMARY_H
#define FASTCAT_FC_RUN_SUMMARY_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_failures.h"
#include "core/fc_read.h"
#include "core/fc_stats.h"
#include "core/fc_stats_engine.h"

typedef struct fc_run_summary fc_run_summary;

typedef struct {
    char *sample_id;
    char *input_type;
    char *barcode_id;
    char *barcode_alias;
    char *run_id;
    char *basecaller_model;
    char *modbase_models;
    char *flow_cell_id;
    char *device_id;
    char *run_date;
    double run_duration_s;
    size_t read_count;
    size_t base_count;
    size_t n5_bp;
    size_t n10_bp;
    size_t n25_bp;
    size_t n50_bp;
    size_t n75_bp;
    size_t n90_bp;
    size_t n95_bp;
    size_t l5;
    size_t l10;
    size_t l25;
    size_t l50;
    size_t l75;
    size_t l90;
    size_t l95;
    uint64_t records_seen;
    uint64_t records_passed_filters;
    uint64_t records_too_long;
    uint64_t records_too_short;
    uint64_t records_low_quality;
    uint64_t records_dust_masked;
} fc_run_summary_row;

fc_run_summary *fc_run_summary_create(size_t barcode_count);
void fc_run_summary_destroy(fc_run_summary *model);

void fc_run_summary_add_read(fc_run_summary *model, size_t barcode, const fc_read *read);
void fc_run_summary_add_filter_event(fc_run_summary *model, size_t barcode, failure_code code);

bool fc_run_summary_seen(const fc_run_summary *model, size_t barcode);

/* True when canonical reads for this barcode report more than one flowcell. */
bool fc_run_summary_has_multiple_flow_cells(const fc_run_summary *model, size_t barcode);

int fc_run_summary_materialize(
    const fc_run_summary *model,
    size_t barcode,
    const fc_stats_engine *engine,
    bool bam_input,
    const read_stats *length_stats,
    const char *sample_id,
    const char *barcode_id,
    fc_run_summary_row *out);

void fc_run_summary_row_destroy(fc_run_summary_row *row);

#endif
