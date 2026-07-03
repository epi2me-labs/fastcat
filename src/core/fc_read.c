#include "core/fc_read.h"

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_fastqcomments.h"
#include "core/fc_filter.h"
#include "core/fc_log.h"

// these are tucked away at the bottom of this file
static void build_tags_string_from_bam(const bam1_t *b, read_meta *meta);
static read_meta *read_meta_create_empty(void);
static void read_meta_destroy(read_meta *meta);
static void read_meta_set_orig_comment(read_meta *meta, const kstring_t *orig_comment);
static void read_meta_reset(read_meta *meta, const kstring_t *orig_comment);
static void read_meta_populate_from_bam(read_meta *meta, const bam1_t *b, const bam_hdr_t *hdr);

fc_read *fc_read_init(void)
{
    fc_read *read = fc_alloc(1, sizeof(*read), "read");
    read->bam = bam_init1();
    if (!read->bam) {
        free(read);
        return NULL;
    }
    return read;
}

void fc_read_destroy(fc_read *read)
{
    if (!read) {
        return;
    }
    if (read->meta) {
        read_meta_destroy(read->meta);
    }
    read->meta = NULL;
    if (read->bam) {
        bam_destroy1(read->bam);
    }
    read->bam = NULL;
    free(read);
}

void fc_read_reset_from_bam(fc_read *read, const bam_hdr_t *hdr, const kstring_t *orig_comment)
{
    if (!read) {
        return;
    }
    if (!read->bam) {
        return;
    }
    if (!read->meta) {
        read->meta = read_meta_create_empty();
    }
    read_meta_reset(read->meta, orig_comment);
    read_meta_populate_from_bam(read->meta, read->bam, hdr);
}

// NOTE: This function is NOT thread-safe due to lazy initialization of tags_str.
// It is designed for single-threaded record processing only.
const kstring_t *fc_read_get_tags_str(const fc_read *read)
{
    if (!read || !read->meta) {
        return NULL;
    }
    if (!read->meta->tags_ready) {
        build_tags_string_from_bam(read->bam, read->meta);
    }
    return read->meta->tags_str;
}

static void read_meta_destroy(read_meta *meta)
{
    if (!meta) {
        return;
    }
    free(meta->orig_comment);
    destroy_rg_info(meta->rg_info);
    free(meta->tags_str->s);
    free(meta->tags_str);
    free(meta);
}

static read_meta *read_meta_create_empty(void)
{
    read_meta *meta = fc_alloc(1, sizeof(read_meta), "meta");
    meta->orig_comment = fc_alloc(1, sizeof(char), "meta->orig_comment");
    meta->orig_comment[0] = '\0';
    meta->tags_str = fc_alloc(1, sizeof(kstring_t), "meta->tags_str");
    ks_initialize(meta->tags_str);
    meta->tags_ready = false;
    return meta;
}

static void read_meta_set_orig_comment(read_meta *meta, const kstring_t *orig_comment)
{
    if (!meta) {
        return;
    }
    size_t len = 0;
    const char *src = NULL;
    if (orig_comment && orig_comment->s && orig_comment->l > 0) {
        len = orig_comment->l;
        src = orig_comment->s;
    }
    // Always realloc to exact size needed (simple and correct)
    size_t need = len + 1;
    meta->orig_comment = fc_realloc(meta->orig_comment, need, "meta->orig_comment");
    if (len > 0) {
        memcpy(meta->orig_comment, src, len);
    }
    meta->orig_comment[len] = '\0';
}

static void read_meta_reset(read_meta *meta, const kstring_t *orig_comment)
{
    read_meta_set_orig_comment(meta, orig_comment);

    destroy_rg_info(meta->rg_info);
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
    meta->mux = 0;
    if (meta->tags_str) {
        ks_clear(meta->tags_str);
    }
    meta->tags_ready = false;
}

static void read_meta_populate_from_bam(read_meta *meta, const bam1_t *b, const bam_hdr_t *hdr)
{
    if (!meta || !b) {
        return;
    }
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
            // Default to unclassified (0) unless we parse a valid barcode number
            meta->ibarcode = 0;
            if (meta->barcode == NULL) {
                // Corrupted aux tag - treat as unclassified rather than crashing
                FC_LOG_WARN("Read '%s' has malformed barcode tag, treating as unclassified.",
                    bam_get_qname(b));
            } else if (strncmp(meta->barcode, "barcode", 7) == 0) {
                const char *num_str = meta->barcode + 7;
                size_t val = 0;
                if (fc_parse_size_strict(num_str, &val) && val < FC_MAX_BARCODES) {
                    meta->ibarcode = val;
                } else {
                    FC_LOG_WARN("Malformed barcode '%s', treating as unclassified", meta->barcode);
                }
            }
        } else if (((strcmp(tag, "BA") == 0) || (strcmp(tag, "AL") == 0)) && type == 'Z') {
            meta->barcode_alias = bam_aux2Z(aux);
        } else if (strcmp(tag, "RN") == 0 && !ierr) {
            meta->read_number = (ival >= 0) ? (size_t)ival : 0;
        } else if (strcmp(tag, "CH") == 0 && !ierr) {
            meta->channel = (ival >= 0) ? (size_t)ival : 0;
        } else if (strcmp(tag, "MX") == 0 && !ierr) {
            meta->mux = (ival >= 0) ? (size_t)ival : 0;
        } else if (strcmp(tag, "NM") == 0 && !ierr) {
            meta->edit_distance = ival;
        } else if ((strcmp(tag, "PI") == 0) && type == 'Z') {
            meta->parent_read = bam_aux2Z(aux);
        } else if (strcmp(tag, "PT") == 0 && !ierr) {
            meta->poly_tail_length = ival;
        } else if (strcmp(tag, "DX") == 0 && !ierr) {
            meta->duplex_code = ival;
        } else if (strcmp(tag, "QS") == 0) {
            errno = 0;
            float qs = (float)bam_aux2f(aux);
            if (!(qs == 0.0f && errno == EINVAL)) {
                meta->mean_q = qs;
            }
        }
        aux = bam_aux_next(b, aux);
    }

    bool good_align = fc_is_primary_qc_pass(b);
    if (good_align && (meta->edit_distance == -1)) {
        FC_LOG_ERROR("Read '%s' does not contain an integer 'NM' tag.", bam_get_qname(b));
        FC_LOG_ERROR("The NM tag is required for all aligned reads (primary, non-supplementary, non-secondary, QC-pass).");
        FC_LOG_ERROR("To fix this:");
        FC_LOG_ERROR("  - Regenerate alignment with a tool that produces NM tags (e.g., minimap2 with default settings)");
        FC_LOG_ERROR("  - Or use 'samtools calmd' to calculate and add NM tags to existing alignments");
        kstring_t rec = { 0, 0, NULL };
        if (sam_format1(hdr, b, &rec) < 0) {
            FC_LOG_ERROR("Failed to format record for error message.");
        } else {
            FC_LOG_ERROR("Offending record: %s", rec.s);
        }
        ks_free(&rec);
        exit(EXIT_FAILURE);
    }

    if (meta->rg && strlen(meta->rg)) {
        meta->rg_info = create_rg_info(meta->rg);
        if (meta->rg_info->runid) {
            meta->runid = meta->rg_info->runid;
        }
        if (meta->rg_info->basecaller) {
            meta->basecaller = meta->rg_info->basecaller;
        }
    }
}

int kseq_to_bam(bam1_t *b, kseq_t *seq)
{
    if (!b || !seq) {
        return -1;
    }

    if (bam_set1(
            b,
            seq->name.l + 1, seq->name.s,
            BAM_FUNMAP, -1, -1, 0,
            0, NULL,
            -1, -1, 0,
            seq->seq.l, seq->seq.s, seq->qual.s,
            0)
        < 0) {
        return -1;
    }

    // bam_set1 uses raw numeric qualities; convert from ASCII
    uint8_t *qual = bam_get_qual(b);
    for (int i = 0; i < b->core.l_qseq; ++i) {
        qual[i] -= 33;
    }

    if (seq->comment.l > 0) {
        if (fastq_comment_to_bam(&seq->comment, b) < 0) {
            return -1;
        }
    }

    return 0;
}

static void build_tags_string_from_bam(const bam1_t *b, read_meta *meta)
{
    if (!meta->tags_str) {
        meta->tags_str = fc_alloc(1, sizeof(kstring_t), "meta->tags_str");
        ks_initialize(meta->tags_str);
    } else {
        ks_clear(meta->tags_str);
    }

    if (!b || sam_format_aux(b, meta->tags_str) < 0) {
        ks_clear(meta->tags_str);
    }
    meta->tags_ready = true;
}
