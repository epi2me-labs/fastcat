#ifndef FASTCAT_FC_SINKS_H
#define FASTCAT_FC_SINKS_H

#include <stdbool.h>
#include <stdio.h>

#include "core/fc_coverage.h"
#include "core/fc_length_qscore_2d.h"
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
    bool use_barcodes;
    size_t barcode_count;
    fc_length_qscore_2d *model;
} fc_length_qscore_2d_sink_opts;

fc_sink fc_length_qscore_2d_sink_create(const fc_length_qscore_2d_sink_opts *opts);

typedef struct {
    const char *output_root;
    bool use_barcodes;
    size_t barcode_count;
    fc_time_bin *model;
} fc_time_bin_sink_opts;

fc_sink fc_time_bin_sink_create(const fc_time_bin_sink_opts *opts);

#endif
