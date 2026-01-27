#ifndef FASTCAT_FC_BAM_UTILS_H
#define FASTCAT_FC_BAM_UTILS_H

#include <stdbool.h>
#include <stdint.h>
#include "htslib/sam.h"

// parameters for bam iteration
typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
    hts_itr_t *iter;
} mplp_data;

/** Set up a bam file for reading records.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *
 *  The return value can be freed with fc_bam_iter_data_destroy.
 *
 */
mplp_data *fc_bam_iter_data_create(
    htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
    const char *chr, hts_pos_t start, hts_pos_t end);

/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void fc_bam_iter_data_destroy(mplp_data *data);

/** Read a bam record (unfiltered).
 *
 *  @param data an mplp_data encoding the bam file to read with filter options.
 *  @param b output pointer.
 *
 */
int fc_bam_read(void *data, bam1_t *b);

/** Create an map of query position to reference position
 *
 *  @param b alignment record
 *
 *  The length of the returned array is b->core->l_qlen.
 */
int *fc_bam_qpos2rpos(bam1_t *b);

uint8_t *fc_bam_get_tag_caseinsensitive(const bam1_t *b, const char tag_name[2]);

#endif
