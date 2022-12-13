#ifndef _BAMSTATS_STATS_H
#define _BAMSTATS_STATS_H

#include <stdbool.h>
#include "htslib/sam.h"

#include "args.h"


/** Generates alignment stats from a region of a bam.
 *
 *  @param fps htsFile pointer
 *  @param idxs hts_idx_t pointer
 *  @param hdrs sam_hdr_t pointer
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param overlap_start whether reads overhanging start should be included.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @param flag_counts size_t flag_counts[8] for output.
 *  @returns void. Prints output to stdout.
 *
 */
void process_bams(
    htsFile *fps, hts_idx_t *idxs, sam_hdr_t *hdrs,
    const char *chr, int start, int end, bool overlap_start,
    const char *read_group, const char tag_name[2], const int tag_value,
    size_t *flag_counts);

#endif
