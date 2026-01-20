#define _GNU_SOURCE
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "thread_pool_internal.h"

#include "../common.h"
#include "../bamcoverage/coverage.h"
#include "../kh_counter.h"
#include "readstats.h"
#include "args.h"
#include "../fc_read.h"
#include "../fc_bam_utils.h"
#include "../fc_filter.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table


// counting alignment flags
// total, primary, ..., unused
static const size_t FLAG_MASK[8] = {
    0, 0, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP, BAM_FQCFAIL, BAM_FDUP, 0
};


/** Initialise flagstat counts struct for a BAM file.
 *
 *  @param n_refs number of reference sequences.
 *
 */
flag_stats* create_flag_stats(size_t n_refs, bool store_unmapped) {
    flag_stats* stats = xalloc(1, sizeof(flag_stats), "flagstat");
    stats->n_refs = n_refs;
    stats->counts = xalloc(n_refs, sizeof(size_t*), "flagstat");
    stats->unmapped = store_unmapped ? xalloc(8, sizeof(size_t), "flagstat") : NULL;

    for (size_t i = 0; i < n_refs; i++) {
        stats->counts[i] = xalloc(8, sizeof(size_t), "flagstat");
    }

    return stats;
}

/** Clean up flagstat counts.
 *
 *  @param stats flagstat counts structure to clean.
 *
 */
void destroy_flag_stats(flag_stats* stats) {
    for (size_t i = 0; i < stats->n_refs; i++) {
        free(stats->counts[i]);
    }
    free(stats->counts);
    free(stats->unmapped);
    free(stats);
}


static inline void process_flagstat_counts(const bam1_t* b, size_t* counts, const int duplex_code ) {
    counts[0] += 1;
    counts[1] += fc_is_primary(b);
    for (size_t i=2; i<6; ++i){
        counts[i] += ((b->core.flag & FLAG_MASK[i]) != 0);
    }
    counts[7] += (duplex_code == 1);
    counts[8] += (duplex_code == -1);
}
// Do all-the-things
void process_bams(
        htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr, const char *sample,
        const char *chr, hts_pos_t start, hts_pos_t end,
        const char *read_group, const char tag_name[2], const int tag_value,
        flag_stats *flag_counts, bool unmapped,
        fc_stats_engine* stats_engine,
        kh_counter_t* runids, kh_counter_t* basecallers, bool force_recalc_qual, cov_writer coverage) {
    if (chr != NULL) {
        if (strcmp(chr, "*") == 0) {
            fprintf(stderr, "Processing: Unplaced reads\n");
        } else {
            fprintf(stderr, "Processing: %s:%zu-%zu\n", chr, (size_t)start, (size_t)end);
        }
    }

    // setup bam reading - reuse our pileup structure, but actually just need iterator
    mplp_data* bam = create_bam_iter_data(
        fp, idx, hdr,
        chr, start, end);
    if (bam == NULL) return;

    int res;
    bam1_t *b = bam_init1();
    char *runid = NULL;
    char *basecaller = NULL;
    char *start_time = NULL;

    while ((res = read_bam(bam, b)) >= 0) {
        if (!fc_filter_int_tag(b, tag_name, tag_value, false)) {
            continue;
        }
        if (!fc_filter_rg(b, read_group, false)) {
            continue;
        }
        // despatch read to coverage calculations
        // NOTE: the writer has its own filters on reads to accept (default 1796: excludes unmapped, secondary, dup, qcfail)
        if (coverage != NULL) {
            coverage_process(coverage, b);
        }

        fc_read* fr = fc_read_from_bam(b, hdr, force_recalc_qual);
        runid = fr->meta->runid ? fr->meta->runid : "";
        basecaller = fr->meta->basecaller ? fr->meta->basecaller : "";
        start_time = fr->meta->start_time ? fr->meta->start_time : "";
        kh_counter_increment(runids, runid);
        kh_counter_increment(basecallers, basecaller);

        char* qname = bam_get_qname(b);
        fc_stats_aligned metrics = {0};
        fc_stats_read_kind kind = fc_stats_engine_add_read(stats_engine, fr, hdr, &metrics);
        if (kind == FC_STATS_UNMAPPED) {
            if (unmapped) {
                // an unmapped read can still have a RNAME and POS, but we
                // ignore that here, because its not a thing we care about
                uint32_t read_length = b->core.l_qseq;
                float mean_quality = fr->meta->mean_q;
                const char *sample_field = sample ? sample : "";
                const char *sample_tab = sample ? "\t" : "";
                fprintf(stdout,
                    "%s\t%s\t%s%s*\tnan\tnan\t" \
                    "nan\tnan\tnan\tnan\t" \
                    "0\t*\t0\t" \
                    "%u\t%.2f\t%s\t" \
                    "0\t0\t0\t0\tnan\tnan\t%d\n",
                    qname, runid, sample_field, sample_tab, //chr, coverage, ref_cover,
                    //qstart, qend, rstart, rend,
                    //aligned_ref_len, direction, length,
                    read_length, mean_quality, start_time,
                    //match, ins, delt, sub, iden, acc
                    fr->meta->duplex_code
                );
                // add to flagstat counts if required
                if (flag_counts != NULL) {
                    process_flagstat_counts(b, flag_counts->unmapped, fr->meta->duplex_code);
                }
            }
            goto FINISH_READ;
        }

        if (flag_counts != NULL) {
            // when we have a target region (as opposed to looping over the whole file),
            // `flag_counts` will only contain one (dynamic) array of counts; otherwise
            // there will be as many dynamic arrays as references in the BAM header
            size_t* counts = (chr != NULL) ? flag_counts->counts[0]
                                           : flag_counts->counts[b->core.tid];
            process_flagstat_counts(b, counts, fr->meta->duplex_code);
        }

        if (kind == FC_STATS_SKIP) {
            goto FINISH_READ;
        }

        {
            const char *ref = (chr != NULL) ? chr : sam_hdr_tid2name(hdr, b->core.tid);
            const char *sample_field = sample ? sample : "";
            const char *sample_tab = sample ? "\t" : "";
            fprintf(stdout,
                "%s\t%s\t%s%s%s\t" \
                "%.4f\t%.4f\t" \
                "%lu\t%lu\t%lu\t%lu\t" \
                "%lu\t%c\t%lu\t%u\t%.2f\t%s\t" \
                "%lu\t%lu\t%lu\t%lu\t%.2f\t%.2f\t%d\n",
                qname, runid, sample_field, sample_tab, ref,
                metrics.coverage, metrics.ref_cover,
                metrics.qstart, metrics.qend, metrics.rstart, metrics.rend,
                metrics.aligned_ref_len, metrics.direction, metrics.length,
                metrics.read_length, metrics.mean_quality, start_time,
                metrics.match, metrics.ins, metrics.delt, metrics.sub,
                metrics.iden, metrics.acc, fr->meta->duplex_code);
        }

FINISH_READ:
        runid = NULL;
        basecaller = NULL;
        start_time = NULL;
        fc_read_destroy(fr);
    }

    destroy_bam_iter_data(bam);
    bam_destroy1(b);

    return;
}
