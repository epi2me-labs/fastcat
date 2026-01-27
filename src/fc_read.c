#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "fc_read.h"
#include "fastqcomments.h"
#include "common.h"
#include "fc_filter.h"

// these two are tucked away at the bottom of this file
static bam1_t* kseq_to_bam(kseq_t* seq);
static void build_tags_string_from_bam(const bam1_t *b, read_meta meta);


fc_read* fc_read_from_kseq(kseq_t* seq) {
    fc_read* r = xalloc(1, sizeof(fc_read), "fc_read");
    bam1_t* bam = kseq_to_bam(seq);
    read_meta meta = fc_read_meta_create(bam, NULL, &seq->comment);

    r->bam = bam;
    r->meta = meta;
    // TODO: honour qs tag if present in the FASTQ comment; currently we match
    // the original fastcat behaviour and always recompute.
    meta->mean_q = mean_qual_naive(seq->qual.s, seq->qual.l);
    r->owns_bam = true;
    return r;
}


fc_read* fc_read_from_bam(bam1_t* bam, const bam_hdr_t* hdr, bool force_recalc_qual) {
    fc_read* r = xalloc(1, sizeof(fc_read), "fc_read");
    r->bam = bam;
    r->owns_bam = false;

    read_meta meta = fc_read_meta_create(bam, hdr, NULL);
    meta->mean_q = resolve_mean_quality(bam, meta, force_recalc_qual);

    r->meta = meta;
    return r;
}


void fc_read_destroy(fc_read* r) {
    if (!r) return;
    if (r->owns_bam && r->bam) {
        bam_destroy1(r->bam);
    }
    if (r->meta) fc_read_meta_destroy(r->meta);
    free(r);
}


read_meta fc_read_meta_create(const bam1_t *b, const bam_hdr_t *hdr, const kstring_t *orig_comment) {
    size_t len = 0;
    const char* src = NULL;
    if (orig_comment && orig_comment->s && orig_comment->l > 0) {
        len = orig_comment->l;
        src = orig_comment->s;
    }

    read_meta meta = xalloc(1, sizeof(_read_meta), "meta");
    meta->orig_comment = xalloc(len + 1, sizeof(char), "meta->orig_comment");
    if (len > 0) {
        memcpy(meta->orig_comment, src, len);
    }
    meta->orig_comment[len] = '\0';
    meta->rg = "";
    meta->rg_info = NULL;
    meta->runid = "";
    meta->basecaller = "";
    meta->flow_cell_id = "";
    meta->barcode = "";
    meta->ibarcode = 0;
    meta->barcode_alias = "";
    meta->start_time = "";
    meta->edit_distance = -1;
    meta->parent_read = "";
    meta->poly_tail_length = -1;
    meta->duplex_code = 0;
    meta->mean_q = -1.0f;
    meta->read_number = 0;
    meta->channel = 0;
    meta->tags_str = xalloc(1, sizeof(kstring_t), "meta->tags_str");
    ks_initialize(meta->tags_str);

    uint8_t *aux = bam_aux_first(b);
    while (aux != NULL) {
        const char *t = bam_aux_tag(aux);
        char tag[3];
        tag[2] = '\0';
        for (int i = 0; i < 2; i++) {
            tag[i] = toupper(t[i]);
        }
        uint8_t type = bam_aux_type(aux);

        int ival = -1;
        bool ierr = false;
        if (type == 'i' || type == 'I' || type == 'c' || type == 'C' || type == 's' || type == 'S') {
            errno = 0;
            ival = bam_aux2i(aux);
            ierr = (ival == 0 && errno == EINVAL);
        }

        if ((strcmp(tag, "RG") == 0) && type == 'Z') {
            meta->rg = bam_aux2Z(aux);
        } else if ((strcmp(tag, "RD") == 0) && type == 'Z') {
            meta->runid = bam_aux2Z(aux);
        } else if ((strcmp(tag, "ST") == 0) && type == 'Z') {
            meta->start_time = bam_aux2Z(aux);
        } else if (((strcmp(tag, "PU") == 0) || (strcmp(tag, "FC") == 0)) && type == 'Z') {
            meta->flow_cell_id = bam_aux2Z(aux);
        } else if (((strcmp(tag, "SM") == 0) || (strcmp(tag, "BC") == 0)) && type == 'Z') {
            meta->barcode = bam_aux2Z(aux);
            if (strcmp(meta->barcode, "unclassified") == 0) {
                meta->ibarcode = 0;
            } else if (strncmp(meta->barcode, "barcode", 7) == 0) {
                meta->ibarcode = (size_t)atoi(meta->barcode + 7);
            } else {
                meta->ibarcode = 0;
            }
        } else if (((strcmp(tag, "BA") == 0) || (strcmp(tag, "AL") == 0)) && type == 'Z') {
            meta->barcode_alias = bam_aux2Z(aux);
        } else if (strcmp(tag, "RN") == 0 && !ierr) {
            meta->read_number = (ival >= 0) ? (size_t)ival : 0;
        } else if (strcmp(tag, "CH") == 0 && !ierr) {
            meta->channel = (ival >= 0) ? (size_t)ival : 0;
        } else if (strcmp(tag, "NM") == 0 && !ierr) {
            meta->edit_distance = ival;
        } else if ((strcmp(tag, "PI") == 0) && type == 'Z') {
            meta->parent_read = bam_aux2Z(aux);
        } else if (strcmp(tag, "PT") == 0 && !ierr) {
            meta->poly_tail_length = ival;
        } else if (strcmp(tag, "DX") == 0 && !ierr) {
            meta->duplex_code = ival;
        } else if (strcmp(tag, "QS") == 0 && (type == 'f' || type == 'd')) {
            errno = 0;
            float qs = bam_aux2f(aux);
            if (!(qs == 0 && errno == EINVAL)) {
                meta->mean_q = qs;
            }
        }
        aux = bam_aux_next(b, aux);
    }

    bool good_align = fc_is_primary_qc_pass(b);
    if (good_align && (meta->edit_distance == -1)) {
        fprintf(stderr, "Read '%s' does not contain an integer 'NM' tag.\n", bam_get_qname(b));
        kstring_t rec = {0, 0, NULL};
        if (sam_format1(hdr, b, &rec) < 0) {
            fprintf(stderr, "Failed to format record for error message.\n");
        } else {
            fprintf(stderr, "%s\n", rec.s);
        }
        ks_free(&rec);
        exit(EXIT_FAILURE);
    }

    if (meta->rg && strlen(meta->rg)) {
        meta->rg_info = create_rg_info(meta->rg);
        if (meta->rg_info->runid) meta->runid = meta->rg_info->runid;
        if (meta->rg_info->basecaller) meta->basecaller = meta->rg_info->basecaller;
    }
    build_tags_string_from_bam(b, meta);
    return meta;
}


void fc_read_meta_destroy(read_meta meta) {
    free(meta->orig_comment);
    destroy_rg_info(meta->rg_info);
    free(meta->tags_str->s);
    free(meta->tags_str);
    free(meta);
}


static bam1_t* kseq_to_bam(kseq_t* seq) {
    bam1_t* b = bam_init1();
    bam_set1(
        b,
        seq->name.l, seq->name.s,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        seq->seq.l, seq->seq.s, seq->qual.s,
        0
    );
    // bam_set1 uses raw numeric qualities; convert from ASCII
    uint8_t* qual = bam_get_qual(b);
    for (int i = 0; i < b->core.l_qseq; ++i) {
        qual[i] -= 33;
    }

    if (seq->comment.l > 0) {
        if (fastq_comment_to_bam(&seq->comment, b) < 0) {
            fprintf(stderr, "Error parsing auxiliary tags for read '%s'\n", seq->name.s);
            exit(1);
        }
    }

    return b;
}


static void build_tags_string_from_bam(const bam1_t *b, read_meta meta) {
    if (!meta->tags_str) {
        meta->tags_str = xalloc(1, sizeof(kstring_t), "meta->tags_str");
        ks_initialize(meta->tags_str);
    } else {
        ks_clear(meta->tags_str);
    }

    if (!b) return;
    if (sam_format_aux(b, meta->tags_str) < 0) {
        ks_clear(meta->tags_str);
    }
}
