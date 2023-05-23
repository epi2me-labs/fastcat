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
#include "bamiter.h"
#include "readstats.h"
#include "args.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table


static const int NOTPRIMARY = BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP;
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
    for (size_t i=b->core.n_cigar - 1; i >= 0; --i){
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


static inline void process_flagstat_counts(bam1_t* b, size_t* counts) {
    counts[0] += 1;
    counts[1] += ((b->core.flag & (NOTPRIMARY)) == 0);
    for (size_t i=2; i<6; ++i){
        counts[i] += ((b->core.flag & FLAG_MASK[i]) != 0);
    }
}


/** Generates alignment stats from a region of a bam.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param sample sample name.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param overlap_start whether reads overhanging start should be included.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @param flag_counts flag_stats pointer.
 *  @param unmapped bool include unmapped reads in output.
 *  @returns void. Prints output to stdout.
 *
 */
void process_bams(
        htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr, const char *sample,
        const char *chr, hts_pos_t start, hts_pos_t end, bool overlap_start,
        const char *read_group, const char tag_name[2], const int tag_value,
        flag_stats *flag_counts, bool unmapped) {
    if (chr != NULL) {
        if (strcmp(chr, "*") == 0) {
            fprintf(stderr, "Processing: Unplaced reads\n");
        } else {
            fprintf(stderr, "Processing: %s:%lld-%lld\n", chr, start, end);
        }
    }

    // setup bam reading - reuse our pileup structure, but actually just need iterator
    mplp_data* bam = create_bam_iter_data(
        fp, idx, hdr,
        chr, start, end, overlap_start,
        read_group, tag_name, tag_value);
    if (bam == NULL) return;

    int res;
    bam1_t *b = bam_init1();
    uint8_t *tag;

    while ((res = read_bam(bam, b) >= 0)) {
        // write a record for unmapped/unplaced
        if (b->core.flag & BAM_FUNMAP){
            if (unmapped) {
                // an unmapped read can still have a RNAME and POS, but we
                // ignore that here, because its not a thing we care about
                char* qname = bam_get_qname(b);
                uint32_t read_length = b->core.l_qseq;
                float mean_quality = mean_qual_from_bam(bam_get_qual(b), read_length);
                if (sample == NULL) {
                    fprintf(stdout,
                        "%s\t*\tnan\tnan\t" \
                        "nan\tnan\tnan\tnan\t" \
                        "0\t*\t0\t%u\t%.3f\t" \
                        "0\t0\t0\t0\tnan\tnan\n",
                        qname, //chr, coverage, ref_cover,
                        //qstart, qend, rstart, rend,
                        //aligned_ref_len, direction, length,
                            read_length, mean_quality
                        //match, ins, delt, sub, iden, acc
                    );
                } else {
                    fprintf(stdout,
                        "%s\t%s\t*\tnan\tnan\t" \
                        "nan\tnan\tnan\tnan\t" \
                        "0\t*\t0\t%u\t%.3f\t" \
                        "0\t0\t0\t0\tnan\tnan\n",
                        qname, sample, //chr, coverage, ref_cover,
                        //qstart, qend, rstart, rend,
                        //aligned_ref_len, direction, length,
                            read_length, mean_quality
                        //match, ins, delt, sub, iden, acc
                    );
                }
                // add to flagstat counts if required
                if (flag_counts != NULL) {
                    process_flagstat_counts(b, flag_counts->unmapped);
                }
            }
            continue;
        }

        if (flag_counts != NULL) {
            // when we have a target region (as opposed to looping over the whole file),
            // `flag_counts` will only contain one (dynamic) array of counts; otherwise
            // there will be as many dynamic arrays as references in the BAM header
            size_t* counts = (chr != NULL) ? flag_counts->counts[0]
                                           : flag_counts->counts[b->core.tid];
            process_flagstat_counts(b, counts);
        }

        // only take "good" primary alignments for further processing
        if (b->core.flag & (NOTPRIMARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        char* qname = bam_get_qname(b);

        // get NM tag
        tag = bam_aux_get((const bam1_t*) b, "NM");
        if (tag == NULL){ // tag isn't present or is currupt
            fprintf(stderr, "Read '%s' does not contain 'NM' tag.\n", qname);
            exit(EXIT_FAILURE);
        }
        int NM = bam_aux2i(tag);
        if (errno == EINVAL) {
            // `get_int_aux_val` returns 0 if setting errno, preventing us from
            // distinguishing between a non-int tag type and an intentional zero. We'll
            // have to check the type ourselves.
            char type = *tag++;
            if (type != 'i' && type != 'I') {
                fprintf(stderr, "Read '%s' contains non-integer 'NM' tag type '%c'\n", qname, type);
                exit(EXIT_FAILURE);
            }
        }
        size_t* stats = create_cigar_stats(b);
        size_t match, ins, delt;
        // some aligners like to get fancy
        match = stats[BAM_CMATCH] + stats[BAM_CEQUAL] + stats[BAM_CDIFF];
        ins = stats[BAM_CINS];
        delt = stats[BAM_CDEL];
        size_t sub = NM - ins - delt;
        size_t length = match + ins + delt;
        float iden = 100 * ((float)(match - sub)) / match;
        float acc = 100 - 100 * ((float)(NM)) / length;
        // we only deal in primary/soft-clipped alignments so length
        // of qseq member is the length of the intact query sequence.
        uint32_t read_length = b->core.l_qseq;
        size_t qstart = get_query_start(b);
        size_t qend = get_query_end(b);
        float mean_quality = mean_qual_from_bam(bam_get_qual(b), read_length);

        float coverage = 100 * ((float)(qend - qstart)) / read_length;
        size_t rstart = b->core.pos;
        size_t rend = bam_endpos(b);
        size_t aligned_ref_len = rend - rstart;
        size_t ref_length = sam_hdr_tid2len(hdr, b->core.tid);
        float ref_cover = 100 * ((float)(aligned_ref_len)) / ref_length;
        char direction = "+-"[bam_is_rev(b)];

        if (sample == NULL) {
            fprintf(stdout,
                "%s\t%s\t" \
                "%.4f\t%.4f\t" \
                "%lu\t%lu\t%lu\t%lu\t" \
                "%lu\t%c\t%lu\t%u\t%.3f\t" \
                "%lu\t%lu\t%lu\t%lu\t%.3f\t%.3f\n",
                qname, (chr != NULL) ? chr : sam_hdr_tid2name(hdr, b->core.tid),
                coverage, ref_cover,
                qstart, qend, rstart, rend,
                aligned_ref_len, direction, length, read_length, mean_quality,
                match, ins, delt, sub, iden, acc);
        } else {
            fprintf(stdout,
                "%s\t%s\t%s\t" \
                "%.4f\t%.4f\t" \
                "%lu\t%lu\t%lu\t%lu\t" \
                "%lu\t%c\t%lu\t%u\t%.3f\t" \
                "%lu\t%lu\t%lu\t%lu\t%.3f\t%.3f\n",
                qname, sample, (chr != NULL) ? chr : sam_hdr_tid2name(hdr, b->core.tid),
                coverage, ref_cover,
                qstart, qend, rstart, rend,
                aligned_ref_len, direction, length, read_length, mean_quality,
                match, ins, delt, sub, iden, acc);
        }
		free(stats);
    }

    destroy_bam_iter_data(bam);
    bam_destroy1(b);

    return;
}
