#ifndef _MODBAMBED_COUNTS_H
#define _MODBAMBED_COUNTS_H

#include <stdbool.h>
#include "htslib/sam.h"

#include "args.h"


/** Generates base counts from a region of a bam.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param overlap_start whether reads overhanging start should be included.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name
 *  @param lowthreshold highest probability to call base as canonical.
 *  @param highthreshold lowest probablity to call base as modified.
 *  @param mod_base BAM code for modified base to report. (e.g. h for 5hmC).
 *  @returns void. Prints output to stdout.
 *
 *
 */
void process_bams(
    htsFile *fps, hts_idx_t *idxs, sam_hdr_t *hdrs,
    const char *chr, int start, int end, bool overlap_start,
    const char *read_group, const char tag_name[2], const int tag_value);


/* Process and print a single region using a threadpool
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param args program arguments.
 *  @param chr reference sequence to process.
 *  @param start reference coordinate to process (0-based).
 *  @param end reference coordiate to process (exclusive).
 *  @param ref reference sequence.
 *
 */
void process_region(
    htsFile **fp, hts_idx_t **idx, sam_hdr_t **hdr,
    arguments_t args, const char *chr, int start, int end, char *ref);

#endif
