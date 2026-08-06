#ifndef FASTCAT_FC_ITER_H
#define FASTCAT_FC_ITER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "core/fc_bam_rg_cache.h"
#include "core/fc_bam_utils.h"
#include "core/fc_common.h"
#include "core/fc_failures.h"
#include "core/fc_filter.h"
#include "core/fc_input_list.h"
#include "core/fc_read.h"
#include "core/fc_regiter.h"

typedef struct fc_run_summary fc_run_summary;

typedef enum {
    FC_ITER_RECORD = 0,
    FC_ITER_EOF = -1,
    FC_ITER_ERROR = -2
} fc_iter_next_status;

typedef struct {
    size_t n_seqs;
    size_t n_bases;
    size_t min_length;
    size_t max_length;
    double mean_quality;
    uint64_t failures[NUM_FAILURE_CODES];
} fc_file_stats;

typedef struct {
    char *buf;
    size_t size;
} dust_buffer;

typedef struct {
    /* common filters */
    size_t min_length;
    size_t max_length;
    float min_qscore;
    bool dust;
    int dust_t;
    int dust_w;
    float max_dust;
    const char *read_group;
    char tag_name[2];
    int tag_value;
    /* BAM-only */
    bool force_recalc_qual;
    const char *region; /* "chr:start-end" */
    const char *bed; /* BED filepath */
    htsThreadPool *thread_pool; /* optional: attach thread pool to BAM handle */
    bool use_barcodes; /* route per-read filter counts by barcode */
    fc_run_summary *run_summary; /* per-barcode filtering counters */
} fc_iter_opts;

/** Disable all read-selection filters while retaining non-filter iterator settings. */
void fc_iter_opts_disable_filters(fc_iter_opts *opts);

typedef struct {
    const char *fname;
    gzFile fp;
    kseq_t *seq;
    fc_iter_opts opts;
    bool apply_filters;
    fc_file_stats stats;
    bool truncated;
    int status;
    double q_sum;
    double q_comp;
    dust_buffer dust;
} fc_fastq_iter;

typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;
    fc_bam_rg_cache *rg_cache;
    mplp_data *bam;
    hts_idx_t *idx;
    fc_iter_opts opts;
    bool apply_filters;
    fc_file_stats stats;
    double q_sum;
    double q_comp;
    bool use_regions;
    regiter reg;
    bool reg_active;
    bool own_idx;
    bool hard_error;
    fc_bam_coordinate_order coordinate_order;
    int status;
    dust_buffer dust;
} fc_bam_iter;

typedef struct {
    char **files;
    size_t n_files;
    size_t idx;
    fc_file_kind expected_kind;
    fc_fastq_iter fqi;
    fc_bam_iter bmi;
    const char *current;
    bool cursor_open;
    bool verbose;
} fc_iter;

void fc_file_stats_init(fc_file_stats *stats);
int fc_iter_init(
    fc_iter *it,
    fc_resolved_inputs *inputs,
    bool verbose);
bool fc_iter_next_file(fc_iter *it, const char **out_file);
int fc_iter_open_record_cursor(
    fc_iter *it,
    const fc_iter_opts *opts);
fc_iter_next_status fc_iter_next_record(fc_iter *it, fc_read *out);
int fc_iter_finish_file(fc_iter *it);
void fc_iter_close(fc_iter *it);

void fc_iter_merge_file_stats(const fc_iter *it, fc_file_stats *dst);
const fc_file_stats *fc_iter_stats(const fc_iter *it);
bool fc_iter_has_hard_error(const fc_iter *it);
#endif
