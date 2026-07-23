#ifndef FASTCAT_FC_SINKS_H
#define FASTCAT_FC_SINKS_H

#include <stdbool.h>
#include <stdio.h>

#include "core/fc_channel_summary.h"
#include "core/fc_coverage.h"
#include "core/fc_length_qscore_2d.h"
#include "core/fc_run_summary.h"
#include "core/fc_speed_qscore_2d.h"
#include "core/fc_stats_engine.h"
#include "core/fc_time_bin.h"
#include "core/fc_writer.h"

typedef struct {
    const char *output_root; // NULL => stream to stdout; non-NULL => write files under <output_root>/{unclassified,barcodeXXXX}
    size_t reads_per_file;
    bool write_bam;
    int threads;
    size_t barcode_count; // size of per-barcode buffers
} fc_output_sink_opts;

fc_sink fc_output_sink_create(const fc_output_sink_opts *opts);

typedef struct {
    const char *output_root;
    const char *sample;
    bool include_alignment;
    bool use_barcodes;
    size_t barcode_count;
} fc_perread_sink_opts;

fc_sink fc_perread_sink_create(const fc_perread_sink_opts *opts);
fc_sink fc_perfile_sink_create(const char *output_root, const char *sample);
typedef struct {
    const char *output_root;
    const char *sample;
    bool use_barcodes;
    size_t barcode_count;
} fc_counter_sink_opts;

fc_sink fc_counter_sink_create(const fc_counter_sink_opts *opts);
typedef struct {
    const char *output_root; // top-level output directory
    bool use_barcodes;
    bool unmapped_override;
    bool write_alignment_metrics;
    bool write_polya;
    fc_stats_engine **engines; // array sized by writer->stats_engine_count
    size_t engine_count;
} fc_hist_sink_opts;

fc_sink fc_hist_sink_create(const fc_hist_sink_opts *opts);

typedef struct {
    bool region_mode;
    const char *region_label;
    const char *output_root;
    const char *sample;
    const bam_hdr_t *hdr;
    bool use_barcodes;
    size_t barcode_count;
} fc_flagstats_sink_opts;

typedef enum {
    FC_FLAGSTAT_TOTAL,
    FC_FLAGSTAT_PRIMARY,
    FC_FLAGSTAT_SECONDARY,
    FC_FLAGSTAT_SUPPLEMENTARY,
    FC_FLAGSTAT_UNMAPPED,
    FC_FLAGSTAT_QCFAIL,
    FC_FLAGSTAT_DUPLICATE,
    FC_FLAGSTAT_DUPLEX,
    FC_FLAGSTAT_DUPLEX_FORMING,
    FC_FLAGSTAT_COUNT
} fc_flagstat_column;

static const char *const fc_flagstat_column_names[FC_FLAGSTAT_COUNT] = {
    "total",
    "primary",
    "secondary",
    "supplementary",
    "unmapped",
    "qcfail",
    "duplicate",
    "duplex",
    "duplex_forming"
};

fc_sink fc_flagstats_sink_create(const fc_flagstats_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool per_base;
    bool use_cigar;
    int exclude_flags;
    int include_flags;
    bool by_strand;
    const bam_hdr_t *hdr;
    const htsThreadPool *pool;
    const char **bed_files;
    const char **bed_names;
    size_t n_beds;
    uint32_t *thresholds;
    size_t n_thresholds;
    uint32_t *segments;
    size_t n_segments;
    bool use_barcodes;
    size_t barcode_count;
} fc_coverage_sink_opts;

fc_sink fc_coverage_sink_create(const fc_coverage_sink_opts *opts);

typedef struct {
    const char *output_root;
    const char *sample;
    bool use_barcodes;
    size_t barcode_count;
    bool bam_input;
    bool unmapped_override;
    fc_stats_engine **engines;
    size_t engine_count;
    fc_run_summary *model;
} fc_run_summary_sink_opts;

fc_sink fc_run_summary_sink_create(const fc_run_summary_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool use_barcodes;
    size_t barcode_count;
    fc_length_qscore_2d *model;
} fc_length_qscore_2d_sink_opts;

fc_sink fc_length_qscore_2d_sink_create(const fc_length_qscore_2d_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool use_barcodes;
    size_t barcode_count;
    fc_speed_qscore_2d *model;
} fc_speed_qscore_2d_sink_opts;

fc_sink fc_speed_qscore_2d_sink_create(
    const fc_speed_qscore_2d_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool use_barcodes;
    size_t barcode_count;
    fc_channel_summary *model;
} fc_channel_summary_sink_opts;

fc_sink fc_channel_summary_sink_create(const fc_channel_summary_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool use_barcodes;
    size_t barcode_count;
    fc_time_bin *model;
} fc_time_bin_sink_opts;

fc_sink fc_time_bin_sink_create(const fc_time_bin_sink_opts *opts);

typedef struct {
    const char *output_root;
    const char *sample;
    bool use_barcodes;
    size_t barcode_count;
    bool bam_input;
    bool unmapped_override;
    bool region_mode;
    const char *region_label;
    const bam_hdr_t *hdr;
    fc_stats_engine **engines;
    size_t engine_count;
    fc_run_summary *run_summary;
    fc_channel_summary *channel_summary;
    fc_length_qscore_2d *length_qscore;
    fc_speed_qscore_2d *speed_qscore;
    fc_time_bin *time_bin;
} fc_summary_json_sink_opts;

fc_sink fc_summary_json_sink_create(const fc_summary_json_sink_opts *opts);

#endif
