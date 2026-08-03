#include "core/fc_stats_engine.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_log.h"

static inline size_t *create_cigar_stats(bam1_t *b)
{
    static const size_t NCODES = 16; // Full 4-bit range; invalid ops 10-15 count to unused slots
    size_t *stats = fc_alloc(NCODES, sizeof(size_t), "read stats");
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        stats[op] += len;
    }
    return stats;
}

static inline size_t get_query_start(bam1_t *b)
{
    uint32_t start_offset = 0;
    uint32_t qlen = b->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if ((start_offset != 0) && (start_offset != qlen)) {
                FC_LOG_FATAL("Invalid clipping in cigar string.");
            }
        } else if (op == BAM_CSOFT_CLIP) {
            start_offset += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }
    return start_offset;
}

static inline size_t get_query_end(bam1_t *b)
{
    uint32_t end_offset = b->core.l_qseq;
    uint32_t qlen = end_offset;
    uint32_t *cigar = bam_get_cigar(b);
    for (int i = b->core.n_cigar - 1; i >= 0; --i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if (end_offset != qlen) {
                FC_LOG_FATAL("Invalid clipping in cigar string.");
            }
        } else if (op == BAM_CSOFT_CLIP) {
            uint32_t oplen = bam_cigar_oplen(cigar[i]);
            if (oplen > end_offset) {
                FC_LOG_FATAL("Read '%s' has soft clip length %u exceeding query length %u.",
                    bam_get_qname(b), oplen, end_offset);
            }
            end_offset -= oplen;
        } else {
            break;
        }
    }
    return end_offset;
}

fc_stats_engine *fc_stats_engine_create(
    bool track_alignment,
    bool track_polya,
    float polya_cover,
    float polya_qual,
    bool polya_rev)
{
    fc_stats_engine *engine = fc_alloc(1, sizeof(fc_stats_engine), "stats engine");
    engine->length_stats = create_length_stats();
    engine->qual_stats = create_qual_stats(QUAL_HIST_WIDTH);

    if (track_alignment) {
        engine->acc_stats = create_qual_stats(ACC_HIST_WIDTH);
        engine->cov_stats = create_qual_stats(COV_HIST_WIDTH);
    }

    if (track_polya) {
        engine->polya_stats = create_length_stats();
    }

    engine->polya_cover = polya_cover;
    engine->polya_qual = polya_qual;
    engine->polya_rev = polya_rev;
    engine->track_alignment = track_alignment;
    return engine;
}

void fc_stats_engine_destroy(fc_stats_engine *engine)
{
    if (engine == NULL) {
        return;
    }
    destroy_length_stats(engine->length_stats);
    destroy_qual_stats(engine->qual_stats);
    destroy_qual_stats(engine->acc_stats);
    destroy_qual_stats(engine->cov_stats);
    destroy_length_stats(engine->polya_stats);
    free(engine);
}

static void fc_stats_engine_add_aligned(fc_stats_engine *engine, const fc_read *fr, const sam_hdr_t *hdr, fc_stats_aligned *out)
{
    if (engine == NULL || fr == NULL || hdr == NULL) {
        return;
    }

    bam1_t *b = fr->bam;
    size_t *stats = create_cigar_stats(b);
    size_t match = stats[BAM_CMATCH] + stats[BAM_CEQUAL] + stats[BAM_CDIFF];
    size_t ins = stats[BAM_CINS];
    size_t delt = stats[BAM_CDEL];
    if (fr->meta->edit_distance < 0) {
        FC_LOG_FATAL("Read '%s' has a negative NM edit distance.", bam_get_qname(b));
    }
    size_t edit_distance = (size_t)fr->meta->edit_distance;
    if (edit_distance < (ins + delt)) {
        FC_LOG_FATAL(
            "Read '%s' has inconsistent NM tag (%zu) vs CIGAR indels (%zu ins + %zu del = %zu).",
            bam_get_qname(b), edit_distance, ins, delt, ins + delt);
    }
    size_t sub = edit_distance - ins - delt;
    size_t length = match + ins + delt;
    if (match == 0 || length == 0) {
        FC_LOG_FATAL("Read '%s' has zero aligned match/length; cannot compute alignment metrics.", bam_get_qname(b));
    }
    float iden = 100 * ((float)(match - sub)) / match;
    float acc = 100 - 100 * ((float)edit_distance) / length;
    if (sub > match || !isfinite(iden) || !isfinite(acc)
        || iden < 0.0 || acc < 0.0 || edit_distance > length) {
        FC_LOG_FATAL("Read '%s' has implausible alignment: sub=%zu match=%zu iden=%.2f acc=%.2f edit_distance=%zu length=%zu",
            bam_get_qname(b), sub, match, iden, acc, edit_distance, length);
    }

    uint32_t read_length = b->core.l_qseq;
    size_t qstart = get_query_start(b);
    size_t qend = get_query_end(b);
    if (qend < qstart || read_length == 0) {
        FC_LOG_FATAL("Read '%s' has invalid query span (%zu-%zu) for read length %u.",
            bam_get_qname(b), qstart, qend, read_length);
    }
    float mean_quality = fr->meta->mean_q;
    float coverage = 100 * ((float)(qend - qstart)) / read_length;
    size_t rstart = b->core.pos;
    size_t rend = bam_endpos(b);
    size_t aligned_ref_len = rend - rstart;
    size_t ref_length = sam_hdr_tid2len(hdr, b->core.tid);
    if (ref_length == 0) {
        FC_LOG_FATAL("Read '%s' is aligned to a zero-length reference.", bam_get_qname(b));
    }
    float ref_cover = 100 * ((float)(aligned_ref_len)) / ref_length;
    if (!isfinite(coverage) || !isfinite(ref_cover)) {
        FC_LOG_FATAL("Read '%s' produced non-finite coverage metrics.", bam_get_qname(b));
    }
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
        *out = (fc_stats_aligned) {
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
    fc_stats_engine *engine, const fc_read *fr,
    const sam_hdr_t *hdr, fc_stats_aligned *out)
{
    if (engine == NULL || fr == NULL || fr->bam == NULL) {
        return FC_STATS_SKIP;
    }

    if (!fc_is_canonical_read_record(fr->bam)) {
        return FC_STATS_SKIP;
    }

    // Generic read metrics are canonical-read metrics: primary mapped and
    // primary unmapped records share one length/quality population.
    add_length_count(engine->length_stats, fr->bam->core.l_qseq);
    add_qual_count(engine->qual_stats, fr->meta->mean_q);

    if (fc_is_unmapped(fr->bam)) {
        return FC_STATS_UNMAPPED;
    }

    if (!fc_is_primary_qc_pass(fr->bam)) {
        return FC_STATS_SKIP;
    }

    if (!engine->track_alignment) {
        FC_LOG_FATAL("Stats engine requires alignment tracking for mapped reads.");
    }
    if (hdr == NULL) {
        FC_LOG_FATAL("Stats engine requires a BAM header for aligned reads.");
    }

    fc_stats_engine_add_aligned(engine, fr, hdr, out);
    return FC_STATS_ALIGNED;
}
