#ifndef _READSTATS_H
#define _READSTATS_H

#include <stdbool.h>

#include "args.h"


/** Generates alignment stats from a region of a bam.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
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
plp_data process_bams(
    const char **bam_file, const char *chr, int start, int end,
    const char *read_group, const char tag_name[2], const int tag_value)


/* Process and print a single region using a threadpool
 *
 * @param args program arguments.
 * @param chr reference sequence to process.
 * @param start reference coordinate to process (0-based).
 * @param end reference coordiate to process (exclusive).
 * @param ref reference sequence.
 *
 */
void process_region(arguments_t args, const char *chr, int start, int end, char *ref);

#endif
