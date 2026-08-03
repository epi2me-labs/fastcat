#include "core/fc_bam_utils.h"

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_log.h"

static const char *format_name(enum htsExactFormat fmt)
{
    switch (fmt) {
    case sam:
        return "SAM";
    case bam:
        return "BAM";
    case cram:
        return "CRAM";
    default:
        return "unknown";
    }
}

static bool bam_headers_match(const sam_hdr_t *left, const sam_hdr_t *right)
{
    if (left->n_targets != right->n_targets) {
        return false;
    }
    for (int tid = 0; tid < left->n_targets; ++tid) {
        const char *left_name = sam_hdr_tid2name(left, tid);
        const char *right_name = sam_hdr_tid2name(right, tid);
        if (strcmp(left_name, right_name) != 0
            || sam_hdr_tid2len(left, tid) != sam_hdr_tid2len(right, tid)) {
            return false;
        }
    }
    return true;
}

fc_bam_input_status fc_bam_input_info_add_header(
    fc_bam_input_info *info, const sam_hdr_t *header)
{
    if (info == NULL || header == NULL || header->n_targets == 0) {
        return FC_BAM_INPUT_OK;
    }
    if (info->reference_header != NULL && !bam_headers_match(info->reference_header, header)) {
        return FC_BAM_INPUT_REFERENCE_MISMATCH;
    }
    if (info->reference_header == NULL) {
        info->reference_header = sam_hdr_dup((sam_hdr_t *)header);
        if (info->reference_header == NULL) {
            return FC_BAM_INPUT_HEADER_COPY_FAILED;
        }
    }
    info->has_references = true;
    return FC_BAM_INPUT_OK;
}

fc_bam_input_status fc_bam_input_preflight(
    const char *const *paths, size_t n_paths, fc_bam_input_info *info,
    const char **failed_path)
{
    if (failed_path != NULL) {
        *failed_path = NULL;
    }
    for (size_t i = 0; i < n_paths; ++i) {
        const char *path = paths[i];
        if (strcmp(path, "-") == 0) {
            info->has_stdin = true;
            continue;
        }
        htsFile *fp = hts_open(path, "rb");
        if (fp == NULL) {
            if (failed_path != NULL) {
                *failed_path = path;
            }
            return FC_BAM_INPUT_OPEN_FAILED;
        }
        sam_hdr_t *header = sam_hdr_read(fp);
        hts_close(fp);
        if (header == NULL) {
            if (failed_path != NULL) {
                *failed_path = path;
            }
            return FC_BAM_INPUT_HEADER_FAILED;
        }
        fc_bam_input_status status = fc_bam_input_info_add_header(info, header);
        sam_hdr_destroy(header);
        if (status != FC_BAM_INPUT_OK) {
            if (failed_path != NULL) {
                *failed_path = path;
            }
            return status;
        }
    }
    return FC_BAM_INPUT_OK;
}

void fc_bam_input_info_destroy(fc_bam_input_info *info)
{
    if (info == NULL) {
        return;
    }
    sam_hdr_destroy(info->reference_header);
    info->reference_header = NULL;
    info->has_references = false;
    info->has_stdin = false;
}

bool fc_bam_header_is_coordinate_sorted(const sam_hdr_t *header)
{
    if (header == NULL) {
        return false;
    }
    kstring_t value = KS_INITIALIZE;
    bool coordinate_sorted = sam_hdr_find_tag_hd((sam_hdr_t *)header, "SO", &value) >= 0
        && strcmp(value.s, "coordinate") == 0;
    free(value.s);
    return coordinate_sorted;
}

bool fc_bam_coordinate_order_check(fc_bam_coordinate_order *order, const bam1_t *record)
{
    if (order == NULL || record == NULL) {
        return false;
    }
    if ((record->core.flag & BAM_FUNMAP) != 0) {
        order->saw_unmapped = true;
        return true;
    }
    int tid = record->core.tid;
    hts_pos_t pos = record->core.pos;
    if (order->saw_unmapped
        || (order->saw_mapped
            && (tid < order->last_mapped_tid
                || (tid == order->last_mapped_tid && pos < order->last_mapped_pos)))) {
        return false;
    }
    order->saw_mapped = true;
    order->last_mapped_tid = tid;
    order->last_mapped_pos = pos;
    return true;
}

void fc_require_bam_input(htsFile *fp, const char *path, const char *operation)
{
    const htsFormat *fmt = hts_get_format(fp);
    enum htsExactFormat exact = fmt ? fmt->format : unknown_format;
    if (exact == bam) {
        return;
    }

    const char *op = operation ? operation : "operation";
    const char *got = format_name(exact);
    if (exact == sam || exact == cram) {
        FC_LOG_FATAL(
            "%s requires BAM input (BGZF-compressed), but '%s' is %s. "
            "Convert first, e.g.: samtools view -b -o reads.bam %s",
            op, path, got, path);
    }
    FC_LOG_FATAL("%s requires BAM input (BGZF-compressed), but '%s' is %s.", op, path, got);
}

/** Set up a bam file for reading records.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  The return value can be freed with fc_bam_iter_data_destroy.
 *
 */
mplp_data *fc_bam_iter_data_create(
    htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
    const char *chr, hts_pos_t start, hts_pos_t end)
{

    mplp_data *data = fc_alloc(1, sizeof(mplp_data), "pileup init data");

    // find the target index for query below
    if (chr == NULL) { // all reads
        data->iter = NULL;
    } else {
        int mytid;
        if (strcmp(chr, "*") == 0) { // unplaced
            mytid = HTS_IDX_NOCOOR;
            start = 0;
            end = INT64_MAX;
        } else {
            mytid = sam_hdr_name2tid(hdr, chr);
            if (mytid < 0) {
                FC_LOG_ERROR("Failed to find reference sequence '%s' in bam.", chr);
                free(data);
                return NULL;
            }
        }
        data->iter = bam_itr_queryi(idx, mytid, start, end);
    }

    // setup bam iterator
    data->fp = fp;
    data->idx = idx;
    data->hdr = hdr;

    return data;
}

bool fc_bam_qualities_missing(const bam1_t *b)
{
    if (NULL == b || b->core.l_qseq <= 0) {
        return true;
    }
    const uint8_t *qual = bam_get_qual(b);
    return NULL == qual || qual[0] == 0xff;
}

/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void fc_bam_iter_data_destroy(mplp_data *data)
{
    bam_itr_destroy(data->iter);
    free(data);
}

/** Read a bam record (unfiltered).
 *
 *  @param data an mplp_data encoding the bam file to read with filter options.
 *  @param b output pointer.
 *
 */
int fc_bam_read(void *data, bam1_t *b)
{
    mplp_data *aux = (mplp_data *)data;
    return aux->iter ? sam_itr_next(aux->fp, aux->iter, b)
                     : sam_read1(aux->fp, aux->hdr, b);
}

/** Create an map of query position to reference position
 *
 *  @param b alignment record
 *
 *  The length of the returned array is b->core->l_qlen.
 */
int *fc_bam_qpos2rpos(const bam1_t *b)
{
    // we only deal in primary/soft-clipped alignments so length
    // of qseq member is the length of the intact query sequence.
    // TODO: add check for alignment being primary / no hard clipping
    uint32_t qlen = b->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(b);
    int *posmap = fc_alloc(qlen, sizeof(uint32_t), "pos_map");
    for (size_t i = 0; i < qlen; ++i) {
        posmap[i] = -1; // unaligned
    }
    int qpos = 0, rpos = b->core.pos;
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        uint32_t take = bam_cigar_type(op);
        if (((take & 0x1) > 0) & ((take & 0x2) > 0)) {
            // consumes query and ref
            for (size_t j = 0; j < len; ++j, ++qpos, ++rpos) {
                posmap[qpos] = rpos;
            }
        } else if ((take & 0x1) > 0) {
            // consumes query only
            qpos += len;
        } else {
            // consumes ref
            rpos += len;
        }
    }
    return posmap;
}

uint8_t *fc_bam_get_tag_caseinsensitive(const bam1_t *b, const char tag_name[2])
{
    if (!tag_name || tag_name[0] == '\0') {
        return NULL;
    }
    char upper_tag[3];
    char lower_tag[3];
    upper_tag[2] = '\0';
    lower_tag[2] = '\0';
    for (int i = 0; i < 2; i++) {
        upper_tag[i] = toupper((unsigned char)tag_name[i]);
        lower_tag[i] = tolower((unsigned char)tag_name[i]);
    }
    uint8_t *ret = bam_aux_get(b, upper_tag);
    if (ret == NULL) {
        ret = bam_aux_get(b, lower_tag);
    }
    return ret;
}

float fc_bam_ensure_mean_quality(bam1_t *b, bool force_recalc_qual)
{
    if (!b) {
        FC_LOG_FATAL("Internal error: null BAM record in mean-quality resolver.");
    }

    // Prefer lowercase qs if both exist.
    const char *target_tag = "qs";
    uint8_t *qs_tag = bam_aux_get(b, "qs");
    if (qs_tag == NULL) {
        qs_tag = bam_aux_get(b, "QS");
        target_tag = "QS";
    }

    if (!force_recalc_qual && qs_tag != NULL) {
        errno = 0;
        float qs = (float)bam_aux2f(qs_tag);
        if (!(qs == 0.0f && errno == EINVAL)) {
            return qs;
        }
    }

    float mean_q = mean_qual_from_bam_naive(bam_get_qual(b), b->core.l_qseq);
    if (qs_tag != NULL) {
        if (bam_aux_update_float(b, target_tag, mean_q) < 0) {
            FC_LOG_FATAL("Failed to update %s tag with recalculated mean quality.", target_tag);
        }
    } else if (bam_aux_append(b, "qs", 'f', sizeof(mean_q), (const uint8_t *)&mean_q) < 0) {
        FC_LOG_FATAL("Failed to append qs tag with recalculated mean quality.");
    }
    return mean_q;
}
