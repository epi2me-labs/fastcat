#ifndef _BAMSTATS_STATS_H
#define _BAMSTATS_STATS_H

#include <stdbool.h>
#include "htslib/sam.h"

#include "args.h"
#include "../stats.h"
#include "../kh_counter.h"


// struct for flagstat counts
typedef struct {
    size_t n_refs;
    size_t** counts;
    size_t* unmapped;
} flag_stats;

/** Create flagstat counts struct for a BAM file.
 *
 *  @param n_refs number of reference sequences.
 *  @param store_unmapped whether to count unmapped reads.
 *
 */
flag_stats* create_flag_stats(size_t n_refs, bool store_unmapped);

/** Clean up flagstat counts.
 *
 *  @param stats flagstat counts structure to clean.
 *
 */
void destroy_flag_stats(flag_stats* stats);


/** Generates alignment stats from a region of a bam.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param sample sample name.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param overlap_start whether reads overhanging start should be included.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @param flag_counts flag_stats pointer.
 *  @param unmapped bool include unmapped reads in output.
 *  @param length_stats read_stats* for accumulating read length information.
 *  @param qual_stats read_stats* for accumulating read quality information.
 *  @param acc_stats read_stats* for accumulating read alignment accuracy information.
 *  @param cov_stats read_stats* for accumulating read alignment coverage information.
 *  @param length_stats_unmapped read_stats* for accumulating read length information for unmapped reads.
 *  @param qual_stats_unmapped read_stats* for accumulating read quality information for unmapped reads.
 *  @param polya_stats read_stats* for accumulating polyA tail length information.
 *  @param polya_cover minimum reference coverage for polyA tail length to be considered.
 *  @param polya_qual minimum mean quality for polyA tail length to be considered.
 *  @param polya_rev whether to allow reverse alignments for polyA tail length.
 *  @param runids kh_counter_t* for accumulating runid information.
 *  @param basecallers kh_counter_t* for accumulating basecaller information.
 *  @param force_recalc_quality whether to recalculate mean quality from phred scores.
 *  @returns void. Prints output to stdout.
 *
 */
void process_bams(
    htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr, const char *sample,
    const char *chr, hts_pos_t start, hts_pos_t end, bool overlap_start,
    const char *read_group, const char tag_name[2], const int tag_value,
    flag_stats *flag_counts, bool unmapped,
    read_stats* length_stats, read_stats* qual_stats, read_stats* acc_stats, read_stats* cov_stats,
    read_stats* length_stats_unmapped, read_stats* qual_stats_unmapped,
    read_stats* polya_stats, float polya_cover, float polya_qual, bool polya_rev,
    kh_counter_t* runids, kh_counter_t* basecallers,
    bool force_recalc_quality);

#endif
