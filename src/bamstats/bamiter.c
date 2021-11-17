#include <errno.h>
#include <string.h>

#include "bamiter.h"
#include "common.h"

/** Set up a bam file for reading (filtered) records.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *
 *  The return value can be freed with destroy_bam_iter_data.
 *
 */
mplp_data *create_bam_iter_data(
        const char *bam_file, const char *chr, int start, int end,
        const char *read_group, const char tag_name[2], const int tag_value) {

    // open bam etc.
    htsFile *fp = hts_open(bam_file, "rb");
    hts_idx_t *idx = sam_index_load(fp, bam_file);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || idx == 0 || fp == 0) {
        hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
        fprintf(stderr, "Failed to read .bam file '%s'.\n", bam_file);
        return NULL;
    }

    // find the target index for query below
    int mytid = -1;
    for (int i=0; i < hdr->n_targets; ++i) {
        if(!strcmp(hdr->target_name[i], chr)) {
            mytid = i;
            break;
        }
    }
    if (mytid == -1) {
        hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
        fprintf(stderr, "Failed to find reference sequence '%s' in bam '%s'.\n", chr, bam_file);
        return NULL;
    }

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->idx = idx; data->hdr = hdr;
    data->iter = bam_itr_queryi(idx, mytid, start, end);
    memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->min_mapQ = 1; data->read_group = read_group;

    return data;
}

/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void destroy_bam_iter_data(mplp_data *data) {
    bam_itr_destroy(data->iter);
    hts_close(data->fp);
    hts_idx_destroy(data->idx);
    sam_hdr_destroy(data->hdr);
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
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        // filter by mapping quality
        if ((int)b->core.qual < aux->min_mapQ) continue;
        // filter by tag
        if (check_tag) {
            tag = bam_aux_get((const bam1_t*) b, aux->tag_name);
            if (tag == NULL){ // tag isn't present or is currupt
                if (aux->keep_missing) {
                    break;
                } else {
                    continue;
                }
            }
            int tag_value = bam_aux2i(tag);
            if (errno == EINVAL) continue; // tag was not integer
            if (tag_value != aux->tag_value) continue;
        }
        // filter by RG (read group):
        if (have_rg) {
            rg = bam_aux_get((const bam1_t*) b, "RG");
            if (rg == NULL) continue;  // missing
            rg_val = bam_aux2Z(rg);
            if (errno == EINVAL) continue;  // bad parse
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
    // ok qseq member is the length of the intact query sequence.
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
