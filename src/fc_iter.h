#ifndef FASTCAT_FC_ITER_H
#define FASTCAT_FC_ITER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <zlib.h>

#include "common.h"
#include "fastcat/parsing.h"
#include "fc_bam_utils.h"
#include "fc_filter.h"
#include "fc_read.h"
#include "htslib/sam.h"
#include "regiter.h"

typedef enum {
    FC_FILE_FASTQ = 0,
    FC_FILE_BAM
} fc_file_kind;

typedef struct {
    char *path;
    fc_file_kind kind;
} fc_file_entry;

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
} fc_iter_opts;

typedef struct {
    const char *fname;
    gzFile fp;
    kseq_t *seq;
    fc_iter_opts opts;
    fc_file_stats stats;
    bool truncated;
    int status;
    double q_sum;
    double q_comp;
} fc_fastq_iter;

int fc_fastq_iter_open(fc_fastq_iter *iter, const char *fname, const fc_iter_opts *opts);
fc_read *fc_fastq_iter_next(fc_fastq_iter *iter);
int fc_fastq_iter_finish(fc_fastq_iter *iter, fc_file_stats *out_stats);

typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;
    mplp_data *bam;
    bam1_t *b;
    hts_idx_t *idx;
    fc_iter_opts opts;
    bool apply_filters;
    fc_file_stats stats;
    double q_sum;
    double q_comp;
    bool has_opts;
    bool use_regions;
    regiter reg;
    bool reg_active;
    bool own_idx;
    int status;
} fc_bam_iter;

int fc_bam_iter_open(
    fc_bam_iter *iter,
    htsFile *fp,
    hts_idx_t *idx,
    sam_hdr_t *hdr,
    const fc_iter_opts *opts);
fc_read *fc_bam_iter_next(fc_bam_iter *iter, const sam_hdr_t *hdr);
int fc_bam_iter_finish(fc_bam_iter *iter, fc_file_stats *out_stats);
void fc_bam_iter_close(fc_bam_iter *iter);

typedef struct {
    fc_file_entry *files;
    size_t n_files;
    size_t idx;
    int recurse;
    fc_fastq_iter fqi;
    fc_bam_iter bmi;
    fc_file_entry *current;
    bool cursor_open;
    htsFile *bam_fp;
    sam_hdr_t *bam_hdr;
    bool verbose;
} fc_iter;

void fc_file_stats_init(fc_file_stats *stats);
int fc_iter_collect_paths(
    char **paths,
    size_t n_paths,
    int recurse,
    bool verbose,
    fc_file_entry **out_files,
    size_t *out_n_files);
int fc_iter_init(fc_iter *it, char **paths, size_t n_paths, int recurse, bool verbose);
bool fc_iter_next_file(fc_iter *it, const fc_file_entry **out_file);
int fc_iter_open_record_cursor(
    fc_iter *it,
    const fc_iter_opts *opts);
fc_read *fc_iter_next_record(fc_iter *it, const sam_hdr_t *hdr);
int fc_iter_finish_file(fc_iter *it);
void fc_iter_close(fc_iter *it);

void fc_iter_merge_file_stats(const fc_iter *it, fc_file_stats *dst);
const fc_file_stats *fc_iter_stats(const fc_iter *it);

#endif
