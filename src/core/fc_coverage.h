#ifndef _BAMCOVERAGE_STATS_H
#define _BAMCOVERAGE_STATS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "core/fc_regiter.h"

typedef struct fc_coverage_summary_target fc_coverage_summary_target;

typedef struct {
    // name used in output file names
    char *name;
    bool whole_chrom;
    // BED data - could simply list complete genome
    bed_regions bed;
    size_t cur_region; // used to detect regions skipped over
    // summary - min, max, total, positions for each region in a BED, and final total entry
    int64_t stats[4];
    int64_t *dist;
    size_t max_cover;
    FILE *fh_summary;
    // distribution - cumulative coverage distribution for each region in a BED, and final total entry
    FILE *fh_dist;
    // thresholds - sparse, user-defined coverage thresholds
    FILE *fh_thresh;
    uint32_t *thresholds;
    size_t n_thresholds;
    fc_coverage_summary_target *summary_target;
    // per-base coverage files
    bool per_base; // whether to write per-base coverage files
    BGZF *fh_fwd;
    char *fn_fwd;
    BGZF *fh_rev;
    char *fn_rev;
    BGZF *fh_tot;
    char *fn_tot;
} _cov_writer_region;

typedef _cov_writer_region *cov_writer_region;

typedef struct {
    char *chrom;
    size_t start;
    size_t end;
    size_t length;
    size_t bases;
    double mean;
    size_t min;
    size_t max;
    double *threshold_fractions;
} fc_coverage_summary_row;

typedef struct {
    char *chrom;
    size_t depth;
    double fraction;
} fc_coverage_distribution_row;

struct fc_coverage_summary_target {
    char *name;
    uint32_t *thresholds;
    size_t n_thresholds;
    fc_coverage_summary_row *summary_rows;
    size_t n_summary_rows;
    size_t cap_summary_rows;
    fc_coverage_distribution_row *distribution_rows;
    size_t n_distribution_rows;
    size_t cap_distribution_rows;
};

typedef struct {
    fc_coverage_summary_target *targets;
    size_t n_targets;
    size_t cap_targets;
} fc_coverage_summary;

typedef struct {
    fc_coverage_summary **by_barcode;
    size_t barcode_count;
} fc_coverage_summaries;

typedef struct {
    // options
    bool use_cigar; // does deletion count
    int exclude_flags; // alignments to exclude
    int include_flags; // alignments to include (applied after)
    // buffers hold deltas
    uint32_t buf_size;
    int32_t *diff_fwd;
    int32_t *diff_rev;
    // current contig state
    const bam_hdr_t *hdr; // cached header
    int32_t tid; // -1 before first record
    int64_t contig_len; // cached from header
    const char *chrom; // cached pointer to hdr->target_name[tid]
    // region handling
    size_t n_beds;
    cov_writer_region *writers; // separate coverage writers for each BED file (and global)
    fc_coverage_summary *summary;
} _cov_writer;

typedef _cov_writer *cov_writer;

void destroy_coverage_writer_region(cov_writer_region writer);

cov_writer init_coverage_writer(
    const char *out_dir, bool write_traces, bool write_tables, bool use_cigar,
    int exclude_flags, int include_flags, bool by_strand,
    const bam_hdr_t *hdr, const htsThreadPool *pool,
    const char **bed_files, const char **bed_names, size_t n_beds,
    uint32_t *thresholds, size_t n_thresholds,
    uint32_t *segments, size_t n_segments,
    fc_coverage_summary *summary);
void destroy_coverage_writer(cov_writer writer);

void coverage_process(cov_writer writer, const bam1_t *b);

fc_coverage_summaries *fc_coverage_summaries_create(size_t barcode_count);
void fc_coverage_summaries_destroy(fc_coverage_summaries *summaries);
fc_coverage_summary *fc_coverage_summaries_prepare(fc_coverage_summaries *summaries, size_t barcode);
const fc_coverage_summary *fc_coverage_summaries_get(const fc_coverage_summaries *summaries, size_t barcode);

#endif // _BAMCOVERAGE_STATS_H
