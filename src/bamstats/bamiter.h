#ifndef _MODBAMBED_BAMITER_H
#define _MODBAMBED_BAMITER_H

#include <stdint.h>
#include <stdbool.h>
#include "htslib/sam.h"

// parameters for bam iteration
typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
    hts_itr_t *iter;
    int min_start;
    int min_mapQ;
    char tag_name[2];
    int tag_value;
    bool keep_missing;
    const char *read_group;
} mplp_data;

/** Set up a bam file for reading (filtered) records.
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
 *  @param tag_value associated with tag_name.
 *
 *  The return value can be freed with destroy_bam_iter_data.
 *
 */
mplp_data *create_bam_iter_data(
    htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
    const char *chr, hts_pos_t start, hts_pos_t end, bool overlap_start,
    const char *read_group, const char tag_name[2], const int tag_value);

/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void destroy_bam_iter_data(mplp_data *data);

/** Read a bam record.
 *
 *  @param data an mplp_data encoding the bam file to read with filter options.
 *  @param b output pointer.
 *
 */
int read_bam(void *data, bam1_t *b);

/** Create an map of query position to reference position
 *
 *  @param b alignment record
 *
 *  The length of the returned array is b->core->l_qlen.
 */
int *qpos2rpos(bam1_t *b);

/** Fetch a BAM tag with case insensitivity
 *
 *  @param b BAM record
 *  @param tag Tag to fetch via bam_aux_get
 *
 */
uint8_t* bam_get_tag_caseinsensitive(const bam1_t* b, char* tag);

#endif
