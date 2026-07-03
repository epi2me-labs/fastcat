#ifndef FASTCAT_FC_REPORT_H
#define FASTCAT_FC_REPORT_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_length_qscore_2d.h"
#include "core/fc_stats_engine.h"
#include "core/fc_time_bin.h"

typedef struct {
    const char *key;
    const char *value;
} fc_report_config_entry;

typedef struct {
    const char *output_root;
    const char *sample;
    const fc_report_config_entry *run_config;
    size_t run_config_count;
    bool use_barcodes;
    size_t barcode_count;
    bool unmapped_override;
    fc_stats_engine **engines;
    size_t engine_count;
    const fc_length_qscore_2d *length_qscore;
    const fc_time_bin *time_bin;
} fc_report_opts;

int fc_write_html_reports(const fc_report_opts *opts);

#endif
