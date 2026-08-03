#ifndef FASTCAT_FC_BAM_UTILS_H
#define FASTCAT_FC_BAM_UTILS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "htslib/sam.h"

#include "core/fc_bam_rg_cache.h"

// parameters for bam iteration
typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
    hts_itr_t *iter;
} mplp_data;

typedef enum {
    FC_BAM_INPUT_OK,
    FC_BAM_INPUT_OPEN_FAILED,
    FC_BAM_INPUT_HEADER_FAILED,
    FC_BAM_INPUT_REFERENCE_MISMATCH,
    FC_BAM_INPUT_HEADER_COPY_FAILED
} fc_bam_input_status;

typedef struct {
    sam_hdr_t *reference_header;
    bool has_references;
    bool has_stdin;
} fc_bam_input_info;

typedef struct {
    int last_mapped_tid;
    hts_pos_t last_mapped_pos;
    bool saw_mapped;
    bool saw_unmapped;
} fc_bam_coordinate_order;

/** Add a BAM header to an input summary.
 *
 * Reference-bearing headers must have matching @SQ names, lengths, and order.
 * The summary owns a copy of the first reference-bearing header.
 */
fc_bam_input_status fc_bam_input_info_add_header(
    fc_bam_input_info *info, const sam_hdr_t *header);

/** Inspect headers for file-backed BAM inputs.
 *
 * Standard input is recorded in info->has_stdin but cannot be opened during
 * preflight. On error, failed_path identifies the file that could not be read
 * or whose reference dictionary did not match.
 */
fc_bam_input_status fc_bam_input_preflight(
    const char *const *paths, size_t n_paths, fc_bam_input_info *info,
    const char **failed_path);

void fc_bam_input_info_destroy(fc_bam_input_info *info);

bool fc_bam_header_is_coordinate_sorted(const sam_hdr_t *header);

/** Record a BAM record in a coordinate-order check.
 *
 * Returns false when mapped records move backwards, or when a mapped record
 * follows an unmapped record.
 */
bool fc_bam_coordinate_order_check(fc_bam_coordinate_order *order, const bam1_t *record);

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

/**
 * Return true when a BAM record has no stored per-base qualities.
 *
 * In BAM, this is represented by a QUAL array filled with 0xff bytes. This is
 * valid for records such as secondary alignments and is distinct from an
 * unreadable record.
 */
bool fc_bam_qualities_missing(const bam1_t *b);

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
int *fc_bam_qpos2rpos(const bam1_t *b);

uint8_t *fc_bam_get_tag_caseinsensitive(const bam1_t *b, const char tag_name[2]);

void fc_require_bam_input(htsFile *fp, const char *path, const char *operation);

/**
 * Resolve mean quality for a BAM record from qs/QS tag unless force_recalc_qual
 * is set, then ensure the effective value is stored back in aux tags.
 *
 * Tag preference is lowercase "qs" then uppercase "QS".
 * If neither exists, a new "qs" tag is appended.
 */
float fc_bam_ensure_mean_quality(bam1_t *b, bool force_recalc_qual);

#endif
