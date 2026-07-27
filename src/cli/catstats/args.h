#ifndef FASTCAT_FC_CATSTATS_ARGS_H
#define FASTCAT_FC_CATSTATS_ARGS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cli/catstats/program_kind.h"
#include "core/fc_iter.h"
#include "core/fc_log.h"

typedef struct {
    // common output + logging
    char *sample;
    char *output_dir;
    bool demultiplex;
    bool per_read;
    bool html_report;
    bool no_legacy_outputs;
    int threads;
    bool verbose;
    bool force_error;
    fc_log_level log_level;
    char *log_file;
    bool log_mono;
    size_t length_qscore_k;
    size_t length_qscore_qgrid;
    size_t length_qscore_qsmooth;
    double time_bin_size_s;
    size_t time_bin_k;

    // common filtering
    size_t min_length;
    size_t max_length;
    float min_qscore;
    bool dust;
    double max_dust;
    size_t dust_w;
    size_t dust_t;
    char *read_group;
    char tag_name[2];
    int tag_value;
    bool force_recalc_qual;

    // fastcat-specific
    int recurse;
    size_t write_bam;
    bool fofn;
    fc_file_kind input_kind;
    size_t reads_per_file;
    fc_path_list input_paths;

    // bamstats-specific
    const char *bam;
    char *region;
    char *bed;

    // coverage options (bamstats)
    bool coverage;
    const char **coverage_beds;
    size_t n_coverage_beds;
    const char **coverage_names;
    size_t n_coverage_names;
    uint32_t *coverage_thresholds;
    size_t n_coverage_thresholds;
    uint32_t *segments;
    size_t n_segments;

    // poly-a options (bamstats)
    bool poly_a;
    float poly_a_cover;
    float poly_a_qual;
    bool poly_a_rev;
} fc_catstats_args;

fc_catstats_args fc_parse_catstats_args(fc_program_kind program, int argc, char **argv);
void fc_destroy_catstats_args(fc_catstats_args *args);

#endif
