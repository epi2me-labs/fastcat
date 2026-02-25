#ifndef FASTCAT_FC_ITER_H
#define FASTCAT_FC_ITER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <zlib.h>

#include "common.h"
#include "fc_failures.h"
#include "fc_bam_utils.h"
#include "fc_filter.h"
#include "fc_read.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "regiter.h"

typedef enum {
    FC_FILE_FASTQ = 0,
    FC_FILE_BAM
} fc_file_kind;

typedef enum {
    FC_ITER_RECORD = 0,
    FC_ITER_EOF = -1,
    FC_ITER_ERROR = -2
} fc_iter_next_status;

bool fc_path_matches_kind(const char *path, fc_file_kind expected_kind);

typedef struct {
    char **paths;
    size_t n_paths;
} fc_path_list;

typedef struct {
    size_t n_seqs;
    size_t n_bases;
    size_t min_length;
    size_t max_length;
    double mean_quality;
    uint64_t failures[NUM_FAILURE_CODES];
} fc_file_stats;

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
    const char *region;      /* "chr:start-end" */
    const char *bed;         /* BED filepath */
    htsThreadPool *thread_pool; /* optional: attach thread pool to BAM handle */
} fc_iter_opts;

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
} fc_fastq_iter;

typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;
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
    int status;
} fc_bam_iter;

typedef struct {
    char **files;
    size_t n_files;
    size_t idx;
    int recurse;
    fc_file_kind expected_kind;
    fc_fastq_iter fqi;
    fc_bam_iter bmi;
    const char *current;
    bool cursor_open;
    bool verbose;
} fc_iter;

void fc_file_stats_init(fc_file_stats *stats);
int fc_iter_collect_paths(
    const fc_path_list *input_paths,
    int recurse,
    fc_file_kind expected_kind,
    bool stdin_is_path_list,
    fc_path_list *out_paths);
int fc_iter_init(
    fc_iter *it,
    const fc_path_list *input_paths,
    int recurse,
    bool verbose,
    fc_file_kind expected_kind,
    bool stdin_is_path_list);
bool fc_iter_next_file(fc_iter *it, const char **out_file);
int fc_iter_open_record_cursor(
    fc_iter *it,
    const fc_iter_opts *opts);
fc_iter_next_status fc_iter_next_record(fc_iter *it, fc_read *out);
int fc_iter_finish_file(fc_iter *it);
void fc_iter_close(fc_iter *it);

void fc_iter_merge_file_stats(const fc_iter *it, fc_file_stats *dst);
const fc_file_stats *fc_iter_stats(const fc_iter *it);
void fc_path_list_free(fc_path_list *paths);

#endif
