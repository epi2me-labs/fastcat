#ifndef FASTCAT_FC_SINKS_H
#define FASTCAT_FC_SINKS_H

#include <stdbool.h>
#include <stdio.h>

#include "bamcoverage/coverage.h"
#include "fc_stats_engine.h"
#include "fc_writer.h"

// barcode 0 is reserved for "unclassified"
#define MAX_BARCODES 1025

typedef struct {
    const char *output_dir;
    size_t reads_per_file;
    bool reheader;
    bool write_bam;
    int threads;
} fc_output_sink_opts;

fc_sink fc_output_sink_create(const fc_output_sink_opts *opts);

typedef struct {
    const char *path;
    const char *sample;
    bool include_alignment;
    bool output_unmapped;
} fc_perread_sink_opts;

fc_sink fc_perread_sink_create(const fc_perread_sink_opts *opts);
fc_sink fc_perfile_sink_create(const char *path, const char *sample);
fc_sink fc_counter_sink_create(const char *runids_path, const char *basecallers_path, const char *sample);
typedef struct {
    const char *hist_dir;
    const char *demux_dir;
    bool use_barcodes;
    bool create_hist_dir;
    bool write_unmapped;
    bool write_alignment_metrics;
    bool write_polya;
    bool unmapped_is_primary;
    fc_stats_engine **engines;
    size_t engine_count;
} fc_hist_sink_opts;

fc_sink fc_hist_sink_create(const fc_hist_sink_opts *opts);

typedef struct flag_stats flag_stats;

typedef struct {
    flag_stats *stats;
    bool region_mode;
    bool track_unmapped;
    const char *out_path;
    const char *sample;
    const bam_hdr_t *hdr;
} fc_flagstats_sink_opts;

fc_sink fc_flagstats_sink_create(const fc_flagstats_sink_opts *opts);

typedef struct {
    const char *out_dir;
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
} fc_coverage_sink_opts;

fc_sink fc_coverage_sink_create(const fc_coverage_sink_opts *opts);

#endif
