#include "fc_stats_engine.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "fc_log.h"
#include "fc_filter.h"

static inline size_t* create_cigar_stats(bam1_t* b) {
    static const size_t NCODES = 10;
    size_t* stats = xalloc(NCODES, sizeof(size_t*), "read stats");
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        stats[op] += len;
    }
    return stats;
}

static inline size_t get_query_start(bam1_t* b) {
    uint32_t start_offset = 0;
    uint32_t qlen = b->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if ((start_offset != 0) && (start_offset != qlen)) {
                FC_LOG_ERROR("Invalid clipping in cigar string.");
                exit(EXIT_FAILURE);
            }
        } else if (op == BAM_CSOFT_CLIP) {
            start_offset += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }
    return start_offset;
}

static inline size_t get_query_end(bam1_t* b) {
    uint32_t end_offset = b->core.l_qseq;
    uint32_t qlen = end_offset;
    uint32_t *cigar = bam_get_cigar(b);
    for (int i = b->core.n_cigar - 1; i >= 0; --i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if (end_offset != qlen) {
                FC_LOG_ERROR("Invalid clipping in cigar string.");
                exit(EXIT_FAILURE);
            }
        } else if (op == BAM_CSOFT_CLIP) {
            end_offset -= bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }
    return end_offset;
}

fc_stats_engine* fc_stats_engine_create(
        bool track_unmapped,
        bool track_alignment,
        bool track_polya,
        float polya_cover,
        float polya_qual,
        bool polya_rev) {
    fc_stats_engine* engine = xalloc(1, sizeof(fc_stats_engine), "stats engine");
    engine->length_stats = create_length_stats();
    engine->qual_stats = create_qual_stats(QUAL_HIST_WIDTH);

    if (track_alignment) {
        engine->acc_stats = create_qual_stats(ACC_HIST_WIDTH);
        engine->cov_stats = create_qual_stats(COV_HIST_WIDTH);
    }

    if (track_unmapped) {
        engine->length_stats_unmapped = create_length_stats();
        engine->qual_stats_unmapped = create_qual_stats(QUAL_HIST_WIDTH);
    }

    if (track_polya) {
        engine->polya_stats = create_length_stats();
    }

    engine->polya_cover = polya_cover;
    engine->polya_qual = polya_qual;
    engine->polya_rev = polya_rev;
    engine->track_alignment = track_alignment;
    engine->track_unmapped = track_unmapped;
    return engine;
}

void fc_stats_engine_destroy(fc_stats_engine* engine) {
    if (engine == NULL) return;
    destroy_length_stats(engine->length_stats);
    destroy_qual_stats(engine->qual_stats);
    destroy_qual_stats(engine->acc_stats);
    destroy_qual_stats(engine->cov_stats);
    destroy_length_stats(engine->length_stats_unmapped);
    destroy_qual_stats(engine->qual_stats_unmapped);
    destroy_length_stats(engine->polya_stats);
    free(engine);
}

static void fc_stats_engine_add_aligned(fc_stats_engine* engine, const fc_read* fr, const sam_hdr_t* hdr, fc_stats_aligned* out) {
    if (engine == NULL || fr == NULL || hdr == NULL) return;

    bam1_t* b = fr->bam;
    size_t* stats = create_cigar_stats(b);
    size_t match = stats[BAM_CMATCH] + stats[BAM_CEQUAL] + stats[BAM_CDIFF];
    size_t ins = stats[BAM_CINS];
    size_t delt = stats[BAM_CDEL];
    size_t sub = fr->meta->edit_distance - ins - delt;
    size_t length = match + ins + delt;
    float iden = 100 * ((float)(match - sub)) / match;
    float acc = 100 - 100 * ((float)(fr->meta->edit_distance)) / length;
    if (iden < 0.0 || acc < 0.0 || (size_t)fr->meta->edit_distance > length) {
        FC_LOG_ERROR("Read '%s' appears to contain implausible alignment information", bam_get_qname(b));
        exit(EXIT_FAILURE);
    }

    uint32_t read_length = b->core.l_qseq;
    size_t qstart = get_query_start(b);
    size_t qend = get_query_end(b);
    float mean_quality = fr->meta->mean_q;
    float coverage = 100 * ((float)(qend - qstart)) / read_length;
    size_t rstart = b->core.pos;
    size_t rend = bam_endpos(b);
    size_t aligned_ref_len = rend - rstart;
    size_t ref_length = sam_hdr_tid2len(hdr, b->core.tid);
    float ref_cover = 100 * ((float)(aligned_ref_len)) / ref_length;
    char direction = "+-"[bam_is_rev(b)];

    if (engine->acc_stats) {
        add_qual_count(engine->acc_stats, acc);
    }
    if (engine->cov_stats) {
        add_qual_count(engine->cov_stats, coverage);
    }

    if (engine->polya_stats != NULL) {
        int polya_len = -1;
        if ((ref_cover >= engine->polya_cover)
                && (!bam_is_rev(b) || engine->polya_rev)
                && mean_quality >= engine->polya_qual) {
            bool has_parent_read = (fr->meta->parent_read != NULL
                && strlen(fr->meta->parent_read) > 0);
            if (!has_parent_read && fr->meta->poly_tail_length >= 0) {
                polya_len = fr->meta->poly_tail_length;
            }
        }
        if (polya_len >= 0) {
            add_length_count(engine->polya_stats, polya_len);
        }
    }

    if (out != NULL) {
        *out = (fc_stats_aligned){
            .read_length = read_length,
            .mean_quality = mean_quality,
            .qstart = qstart,
            .qend = qend,
            .rstart = rstart,
            .rend = rend,
            .aligned_ref_len = aligned_ref_len,
            .ref_length = ref_length,
            .match = match,
            .ins = ins,
            .delt = delt,
            .sub = sub,
            .length = length,
            .coverage = coverage,
            .ref_cover = ref_cover,
            .iden = iden,
            .acc = acc,
            .direction = direction
        };
    }

    free(stats);
}

fc_stats_read_kind fc_stats_engine_add_read(
        fc_stats_engine* engine, const fc_read* fr,
        const sam_hdr_t* hdr, fc_stats_aligned* out) {
    if (engine == NULL || fr == NULL || fr->bam == NULL) return FC_STATS_SKIP;

    if (fc_is_unmapped(fr->bam)) {
        if (engine->track_unmapped) {
            add_length_count(engine->length_stats_unmapped, fr->bam->core.l_qseq);
            add_qual_count(engine->qual_stats_unmapped, fr->meta->mean_q);
        }
        return FC_STATS_UNMAPPED;
    }

    if (!fc_is_primary_qc_pass(fr->bam)) {
        return FC_STATS_SKIP;
    }

    if (!engine->track_alignment) {
        FC_LOG_ERROR("Stats engine requires alignment tracking for mapped reads.");
        exit(EXIT_FAILURE);
    }
    if (hdr == NULL) {
        FC_LOG_ERROR("Stats engine requires a BAM header for aligned reads.");
        exit(EXIT_FAILURE);
    }

    add_length_count(engine->length_stats, fr->bam->core.l_qseq);
    add_qual_count(engine->qual_stats, fr->meta->mean_q);
    fc_stats_engine_add_aligned(engine, fr, hdr, out);
    return FC_STATS_ALIGNED;
}
