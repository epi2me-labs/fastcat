#include <ctype.h>
#include <errno.h>
#include <string.h>

#include "bamiter.h"
#include "common.h"

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
        const char *read_group, const char tag_name[2], const int tag_value) {

    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");

    // find the target index for query below
    if (chr == NULL) {  // all reads
        data->iter = NULL;
    } else {
        int mytid;
        if (strcmp(chr, "*") == 0) { // unplaced
            mytid = HTS_IDX_NOCOOR;
            start = 0; end = INT64_MAX;
        } else {
            mytid = sam_hdr_name2tid(hdr, chr);
            if (mytid < 0) {
                fprintf(stderr, "Failed to find reference sequence '%s' in bam.\n", chr);
                free(data);
                return NULL;
            }
        }
        data->iter = bam_itr_queryi(idx, mytid, start, end);
    }

    // setup bam interator
    data->fp = fp; data->idx = idx; data->hdr = hdr;
    data->min_start = overlap_start ? -1 : start; // unmapped reads have pos -1
    memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->min_mapQ = 0; data->read_group = read_group;

    return data;
}

/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void destroy_bam_iter_data(mplp_data *data) {
    bam_itr_destroy(data->iter);
    free(data);
}


/** Read a bam record.
 *
 *  @param data an mplp_data encoding the bam file to read with filter options.
 *  @param b output pointer.
 *
 */
int read_bam(void *data, bam1_t *b) {
    mplp_data *aux = (mplp_data*) data;
    uint8_t *tag;
    bool check_tag = (strcmp(aux->tag_name, "") != 0);
    bool have_rg = (aux->read_group != NULL);
    uint8_t *rg;
    char *rg_val;
    int ret;
    while (1) {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if (ret<0) break;
        // only take primary alignments
        //if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        // maybe remove reads overlapping start
        if (b->core.pos < aux->min_start) continue;
        // filter by mapping quality
        if ((int)b->core.qual < aux->min_mapQ) continue;
        // filter by tag
        if (check_tag) {
            tag = bam_get_tag_caseinsensitive((const bam1_t*) b, aux->tag_name);
            if (tag == NULL){ // tag isn't present or is currupt
                if (aux->keep_missing) {
                    break;
                } else {
                    continue;
                }
            }
            int tag_value = bam_aux_tag_int(tag);
            if (tag_value == 0 && errno == EINVAL) continue; // tag was not integer
            if (tag_value != aux->tag_value) continue;
        }
        // filter by RG (read group):
        if (have_rg) {
            rg = bam_get_tag_caseinsensitive((const bam1_t*) b, "RG");
            if (rg == NULL) continue;  // missing
            rg_val = bam_aux2Z(rg);
            if (rg_val == 0 && errno == EINVAL) continue;  // bad parse
            if (strcmp(aux->read_group, rg_val) != 0) continue;  // not wanted
        }
        break;
    }
    return ret;
}


/** Create an map of query position to reference position
 *
 *  @param b alignment record
 *
 *  The length of the returned array is b->core->l_qlen.
 */
int *qpos2rpos(bam1_t *b) {
    // we only deal in primary/soft-clipped alignments so length
    // of qseq member is the length of the intact query sequence.
    // TODO: add check for alignment being primary / no hard clipping
    uint32_t qlen = b->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(b);
    int *posmap = xalloc(qlen, sizeof(uint32_t), "pos_map");
    for (size_t i = 0; i < qlen; ++i) posmap[i] = -1;  // unaligned
    int qpos = 0, rpos = b->core.pos;
    for (size_t i = 0; i < b->core.n_cigar; ++i){
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        uint32_t take = bam_cigar_type(op);
        if (((take&0x1)>0) & ((take&0x2)>0)) {
            // consumes query and ref
            for (size_t j = 0; j < len; ++j, ++qpos, ++rpos) {
                posmap[qpos] = rpos;
            }
        }
        else if ((take&0x1)>0) {
            // consumes query only
            qpos += len;
        }
        else {
            // consumes ref
            rpos += len;
        }
    }
    return posmap;
}

/** Fetch a BAM tag with case insensitivity
 *
 *  @param b BAM record
 *  @param tag Tag to fetch via bam_aux_get
 *
 */
uint8_t* bam_get_tag_caseinsensitive(const bam1_t* b, char* tag) {

    uint8_t* ret;
    char upper_tag[3];
    char lower_tag[3];
    upper_tag[2] = '\0';
    lower_tag[2] = '\0';
    for (int i = 0; i < 2; i++) {
        upper_tag[i] = toupper(tag[i]);
        lower_tag[i] = tolower(tag[i]);
    }
    // Try uppercase variant
    ret = bam_aux_get((const bam1_t*) b, upper_tag);
    if (ret == NULL){
        // Try lowercase variant
        ret = bam_aux_get((const bam1_t*) b, lower_tag);
    }
    return ret;
}

/** Translate an integer typed aux tag to integer form after zeroing errno.
 *
 *  @param tag Aux tag data to translate to int
 *
 */
int bam_aux_tag_int(uint8_t* tag) {
    // get_int_aux_val (inside bam_aux2i) unhelpfully returns zero if the
    // provided tag data is invalid. this is unfortunate, as zero is a valid
    // value for an integer type field.
    // it is not usually proper to manually reset errno, however it is simply
    // not possible to discern between a true zero value and an invalid zero
    // that requires handling of EINVAL. we must reset errno to prevent any
    // previously encountered EINVAL errors triggering behaviour on a true zero.
    errno = 0;
    // it is now safe to check ret == 0 && errno == EINVAL
    return bam_aux2i(tag);
}
