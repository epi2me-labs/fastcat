#include "fc_bam_utils.h"

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#include "common.h"
#include "fc_log.h"

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
