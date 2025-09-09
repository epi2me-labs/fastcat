#ifndef _BAMCOVERAGE_STATS_H
#define _BAMCOVERAGE_STATS_H

#include <stdint.h>
#include <stdbool.h>
#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "regiter.h"

typedef struct {
    // name used in output file names
    char* name;
    bool whole_chrom;
    // BED data - could simply list complete genome
    bed_regions bed;
    size_t cur_region; // used to detect regions skipped over
    // summary - min, max, total, positions for each region in a BED, and final total entry
    int64_t stats[4];
    int64_t* dist;
    size_t max_cover;
    FILE* fh_summary;
    // distribution - cumulative coverage distribution for each region in a BED, and final total entry
    FILE* fh_dist;
    // thresholds - sparse, user-defined coverage thresholds
    FILE* fh_thresh;
    uint32_t* thresholds;
    size_t n_thresholds;
    // per-base coverage files
    bool per_base; // whether to write per-base coverage files
    BGZF* fh_fwd;
    char* fn_fwd;
    BGZF* fh_rev;
    char* fn_rev;
    BGZF* fh_tot;
    char* fn_tot;
} _cov_writer_region;

typedef _cov_writer_region* cov_writer_region;


typedef struct {
    // options
    bool use_cigar;    // does deletion count
    int exclude_flags; // alignments to exclude
    int include_flags; // alignments to include (applied after)
    // buffers hold deltas
    uint32_t buf_size;
    int32_t* diff_fwd;
    int32_t* diff_rev;
    // current contig state
    const bam_hdr_t* hdr;    // cached header
    int32_t  tid;            // -1 before first record
    int64_t  contig_len;     // cached from header
    const char* chrom;       // cached pointer to hdr->target_name[tid]
    // region handling
    size_t n_beds;
    cov_writer_region* writers;  // separate coverage writers for each BED file (and global)
} _cov_writer;

typedef _cov_writer* cov_writer;


cov_writer_region init_coverage_writer_region(
        const char* out_dir, const char* name, bool per_base, bool by_strand, const char* bed_fname,
        uint32_t* thresholds, size_t n_thresholds, uint32_t segment_length,
        const bam_hdr_t* hdr, const htsThreadPool* pool);
void destroy_coverage_writer_region(cov_writer_region writer);


cov_writer init_coverage_writer(
        const char* out_dir, bool per_base, bool use_cigar,
        int exclude_flags, int include_flags, bool by_strand,
        const bam_hdr_t* hdr, const htsThreadPool* pool,
        char** bed_files, char** bed_names, size_t n_beds,
        uint32_t* thresholds, size_t n_thresholds,
        uint32_t* segments, size_t n_segments);
void destroy_coverage_writer(cov_writer writer);


void coverage_process(cov_writer writer, const bam1_t* b);

#endif // _BAMCOVERAGE_STATS_H
