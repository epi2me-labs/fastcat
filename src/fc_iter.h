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

typedef struct {
    size_t n_seqs;
    size_t n_bases;
    size_t min_length;
    size_t max_length;
    double mean_quality;
    uint64_t failures[NUM_FAILURE_CODES];
} fc_file_stats;

typedef struct {
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
} fc_fastq_iter_opts;

typedef struct {
    const char *read_group;
    char tag_name[2];
    int tag_value;
    bool force_recalc_qual;
    size_t min_length;
    size_t max_length;
    float min_qscore;
    bool dust;
    int dust_t;
    int dust_w;
    float max_dust;
} fc_bam_iter_opts;

typedef struct {
    const char *fname;
    gzFile fp;
    kseq_t *seq;
    fc_fastq_iter_opts opts;
    fc_filter_opts filter_opts;
    fc_file_stats stats;
    bool truncated;
    int status;
    size_t n_valid;
    size_t total_bases;
    size_t min_len;
    size_t max_len;
    double average_q;
    double q_sum;
} fc_fastq_iter;

int fc_fastq_iter_open(fc_fastq_iter *iter, const char *fname, const fc_fastq_iter_opts *opts);
fc_read *fc_fastq_iter_next(fc_fastq_iter *iter);
int fc_fastq_iter_finish(fc_fastq_iter *iter, fc_file_stats *out_stats);

typedef struct {
    mplp_data *bam;
    bam1_t *b;
    fc_bam_iter_opts opts;
    fc_filter_opts filter_opts;
    bool apply_filters;
    fc_file_stats stats;
    size_t n_valid;
    size_t total_bases;
    size_t min_len;
    size_t max_len;
    double average_q;
    double q_sum;
    bool has_opts;
    int status;
} fc_bam_iter;

int fc_bam_iter_open(
    fc_bam_iter *iter,
    htsFile *fp,
    hts_idx_t *idx,
    sam_hdr_t *hdr,
    const char *chr,
    hts_pos_t start,
    hts_pos_t end,
    const fc_bam_iter_opts *opts);
fc_read *fc_bam_iter_next(fc_bam_iter *iter, const sam_hdr_t *hdr);
int fc_bam_iter_finish(fc_bam_iter *iter, fc_file_stats *out_stats);
void fc_bam_iter_close(fc_bam_iter *iter);

#endif
