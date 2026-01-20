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
#include "../stats.h"
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

// Count number of each cigar operation in an alignment
inline size_t* create_cigar_stats (bam1_t* b) {
	static const size_t NCODES = 10;
    size_t* stats = xalloc(NCODES, sizeof(size_t*), "read stats");
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i){
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
		stats[op] += len;
	}
	return stats;
}

// Find the first aligned position of the query sequence
inline size_t get_query_start (bam1_t* b) {
    uint32_t start_offset = 0;
    uint32_t qlen = b->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i){
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if ((start_offset != 0) && (start_offset != qlen)) {
                fprintf(stderr, "Invalid clipping in cigar string.\n");
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

// Find the last aligned position of the query sequence
inline size_t get_query_end(bam1_t* b) {
    // TODO: assume l_qseq correct
    uint32_t end_offset = b->core.l_qseq;
    uint32_t qlen = end_offset;
    uint32_t *cigar = bam_get_cigar(b);
    for (int i=b->core.n_cigar - 1; i >= 0; --i){
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            if (end_offset != qlen) {
                fprintf(stderr, "Invalid clipping in cigar string.\n");
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
        read_stats* length_stats, read_stats* qual_stats, read_stats* acc_stats, read_stats* cov_stats,
        read_stats* length_stats_unmapped, read_stats* qual_stats_unmapped,
        read_stats* polya_stats, float polya_cover, float polya_qual, bool polya_rev,
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

        // write a record for unmapped/unplaced
        if (fc_is_unmapped(b)) {
            if (unmapped) {
                // an unmapped read can still have a RNAME and POS, but we
                // ignore that here, because its not a thing we care about
                char* qname = bam_get_qname(b);
                uint32_t read_length = b->core.l_qseq;
                float mean_quality = fr->meta->mean_q;
                if (sample == NULL) {
                    fprintf(stdout,
                        "%s\t%s\t*\tnan\tnan\t" \
                        "nan\tnan\tnan\tnan\t" \
                        "0\t*\t0\t" \
                        "%u\t%.2f\t%s\t" \
                        "0\t0\t0\t0\tnan\tnan\t%d\n",
                        qname, runid, //chr, coverage, ref_cover,
                        //qstart, qend, rstart, rend,
                        //aligned_ref_len, direction, length,
                        read_length, mean_quality, start_time,
                        //match, ins, delt, sub, iden, acc
                        fr->meta->duplex_code
                    );
                } else {
                    fprintf(stdout,
                        "%s\t%s\t%s\t*\tnan\tnan\t" \
                        "nan\tnan\tnan\tnan\t" \
                        "0\t*\t0\t" \
                        "%u\t%.2f\t%s\t" \
                        "0\t0\t0\t0\tnan\tnan\t%d\n",
                        qname, runid, sample, //chr, coverage, ref_cover,
                        //qstart, qend, rstart, rend,
                        //aligned_ref_len, direction, length,
                        read_length, mean_quality, start_time,
                        //match, ins, delt, sub, iden, acc
                        fr->meta->duplex_code
                    );
                }
                // add to flagstat counts if required
                if (flag_counts != NULL) {
                    process_flagstat_counts(b, flag_counts->unmapped, fr->meta->duplex_code);
                }

                // accumulate stats into histogram
                add_length_count(length_stats_unmapped, read_length);
                add_qual_count(qual_stats_unmapped, mean_quality);
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

        // only take "good" primary alignments for further processing
        if (!fc_is_primary_qc_pass(b)) {
            goto FINISH_READ;
        }
        char* qname = bam_get_qname(b);

        size_t* stats = create_cigar_stats(b);
        size_t match, ins, delt;
        // some aligners like to get fancy
        match = stats[BAM_CMATCH] + stats[BAM_CEQUAL] + stats[BAM_CDIFF];
        ins = stats[BAM_CINS];
        delt = stats[BAM_CDEL];
        size_t sub = fr->meta->edit_distance - ins - delt;
        size_t length = match + ins + delt;
        float iden = 100 * ((float)(match - sub)) / match;
        float acc = 100 - 100 * ((float)(fr->meta->edit_distance)) / length;
        // some things we've seen go wrong
        // explode now because there is almost certainly something wrong with the tags
        // and calling add_qual_count with a value less than zero will cause a segfault
        if (iden < 0.0 || acc < 0.0 || (size_t)fr->meta->edit_distance > length) {
            fprintf(stderr, "Read '%s' appears to contain implausible alignment information\n", qname);
            exit(EXIT_FAILURE);
        }
        // we only deal in primary/soft-clipped alignments so length
        // of qseq member is the length of the intact query sequence.
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

        // accumulate stats into histogram
        add_length_count(length_stats, read_length);
        add_qual_count(qual_stats, mean_quality);
        add_qual_count(acc_stats, acc);
        add_qual_count(cov_stats, coverage);

        // get poly-A tail length. For now we require:
        //    i) "good" coverage on reference, i.e. "full length"
        //   ii) read is sense strand, i.e. fwd alignment
        //  iii) "good" mean quality
        //   iv) no split reads
        if (polya_stats != NULL) {
            int polya_len = -1;
            if ((ref_cover >= polya_cover)
                    && (!bam_is_rev(b) || polya_rev)
                    && mean_quality >= polya_qual) {
                bool has_parent_read = (fr->meta->parent_read != NULL
                    && strlen(fr->meta->parent_read) > 0);
                if (!has_parent_read && fr->meta->poly_tail_length >= 0) {
                    polya_len = fr->meta->poly_tail_length;
                }
            }
            if (polya_len >= 0) {
                add_length_count(polya_stats, polya_len);
            }
        }

        if (sample == NULL) {
            fprintf(stdout,
                "%s\t%s\t%s\t" \
                "%.4f\t%.4f\t" \
                "%lu\t%lu\t%lu\t%lu\t" \
                "%lu\t%c\t%lu\t%u\t%.2f\t%s\t" \
                "%lu\t%lu\t%lu\t%lu\t%.2f\t%.2f\t%d\n",
                qname, runid, (chr != NULL) ? chr : sam_hdr_tid2name(hdr, b->core.tid),
                coverage, ref_cover,
                qstart, qend, rstart, rend,
                aligned_ref_len, direction, length, read_length, mean_quality, start_time,
                match, ins, delt, sub, iden, acc, fr->meta->duplex_code);
        } else {
            fprintf(stdout,
                "%s\t%s\t%s\t%s\t" \
                "%.4f\t%.4f\t" \
                "%lu\t%lu\t%lu\t%lu\t" \
                "%lu\t%c\t%lu\t%u\t%.2f\t%s\t" \
                "%lu\t%lu\t%lu\t%lu\t%.2f\t%.2f\t%d\n",
                qname, runid, sample, (chr != NULL) ? chr : sam_hdr_tid2name(hdr, b->core.tid),
                coverage, ref_cover,
                qstart, qend, rstart, rend,
                aligned_ref_len, direction, length, read_length, mean_quality, start_time,
                match, ins, delt, sub, iden, acc, fr->meta->duplex_code);
        }
		free(stats);

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
