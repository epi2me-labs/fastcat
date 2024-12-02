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
#include "../stats.h"
#include "../kh_counter.h"
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
    counts[1] += ((b->core.flag & (NOTPRIMARY)) == 0);
    for (size_t i=2; i<6; ++i){
        counts[i] += ((b->core.flag & FLAG_MASK[i]) != 0);
    }
    counts[7] += (duplex_code == 1);
    counts[8] += (duplex_code == -1);
}


// see section 4.2.4 of the SAM spec for more details
#define IS_INTEGER_TAG(t) ((t) == 'i' || (t) == 'I' || (t) == 'c' || (t) == 'C' || (t) == 's' || (t) == 'S')

#define N_TAGS 8
typedef struct {
    char *RG;  // read group
    char *RD;  // read group (old skool)
    char *st;  // start time
    int NM;    // edit distance
    int pi;    // parent read
    int pt;    // poly-t/a tail length
    float qs;  // quality score
    int dx;    // duplex
} bam_tags_t;


// Function to fetch tags from a bam1_t record
bam_tags_t fetch_bam_tags(const bam1_t *b, const bam_hdr_t *header) {
    // default duplex tag to simple read, everything else as invalid
    bam_tags_t tags = {NULL, NULL, NULL, -1, -1, -1, 0, -1};  

    uint8_t *aux = bam_aux_first(b);
    int n_tags = 0;
    while (n_tags < 100 && aux != NULL) {
        n_tags++;
        const char *t = bam_aux_tag(aux);
        char tag[3];
        tag[2] = '\0';
        for (int i = 0; i < 2; i++) {
            tag[i] = toupper(t[i]);
        }
        uint8_t type = bam_aux_type(aux);

        // do this here to avoid repeating below
        int ival = -1;
        bool ierr = false;
        if (IS_INTEGER_TAG(type)) {
            errno = 0;
            ival = bam_aux2i(aux);
            ierr = (ival == 0 && errno == EINVAL);
        }
        
        if ((strcmp(tag, "RG") == 0) && (tags.RG == NULL) && type == 'Z') {
            tags.RG = strdup(bam_aux2Z(aux));
        } else if ((strcmp(tag, "RD") == 0) && (tags.RD == NULL) && type == 'Z') {
            tags.RD = strdup(bam_aux2Z(aux));
        } else if ((strcmp(tag, "ST") == 0) && (tags.st == NULL) && type == 'Z') {
            tags.st = strdup(bam_aux2Z(aux));
        } else if (strcmp(tag, "NM") == 0 && !ierr) {
            tags.NM = ival;
        } else if (strcmp(tag, "PI") == 0 && !ierr) {
            tags.pi = ival;
        } else if (strcmp(tag, "PT") == 0 && !ierr) {
            tags.pt = ival;
        } else if (strcmp(tag, "DX") == 0 && !ierr) {
            tags.dx = ival;
        } else if (strcmp(tag, "QS") == 0 && (type == 'f')) {
            errno = 0;
            tags.qs = bam_aux2f(aux);
            if (tags.qs == 0 && errno == EINVAL) {
                tags.qs = -1;
            }
        }
        else {
            // we added above when we shouldn't have
            n_tags--;
        }

        aux = bam_aux_next(b, aux);
    }

    // Check we have all the tags we need
    // note theres weird corner case of duplicate tags, but when does that happen?
    bool good_align = ((b->core.flag & (NOTPRIMARY | BAM_FQCFAIL | BAM_FDUP)) == 0);
    if (good_align && (tags.NM == -1)) {
        fprintf(stderr, "Read '%s' does not contain an integer 'NM' tag.\n", bam_get_qname(b));
        kstring_t rec = {0, 0, NULL};
        if (sam_format1(header, b, &rec) < 0) {
            fprintf(stderr, "Failed to format record for error message.\n");
        } else {
            fprintf(stderr, "%s\n", rec.s);
        }
        ks_free(&rec);
        exit(EXIT_FAILURE);
    }
    return tags;
}


void free_bam_tags(bam_tags_t *tags) {
    if (!tags) return;
    if (tags->RG) free(tags->RG);
    if (tags->RD) free(tags->RD);
    if (tags->st) free(tags->st);
    //free(tags);  // we stack allocate all these now
}


// Do all-the-things
void process_bams(
        htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr, const char *sample,
        const char *chr, hts_pos_t start, hts_pos_t end, bool overlap_start,
        const char *read_group, const char tag_name[2], const int tag_value,
        flag_stats *flag_counts, bool unmapped,
        read_stats* length_stats, read_stats* qual_stats, read_stats* acc_stats, read_stats* cov_stats,
        read_stats* length_stats_unmapped, read_stats* qual_stats_unmapped,
        read_stats* polya_stats, float polya_cover, float polya_qual, bool polya_rev,
        kh_counter_t* runids, kh_counter_t* basecallers,
        bool force_recalc_qual) {
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
        chr, start, end, overlap_start,
        read_group, tag_name, tag_value);
    if (bam == NULL) return;

    int res;
    bam1_t *b = bam_init1();
    readgroup* rg_info = NULL;
    char *runid = NULL;
    char *basecaller = NULL;
    char *start_time = NULL;

    while ((res = read_bam(bam, b) >= 0)) {
        // get all our tags
        bam_tags_t tags = fetch_bam_tags(b, hdr);

        // get info from readgroup, note we could use subitems from readgroup
        // here more directly, but this is to be consistent with fastcat where
        // we only have the readgroup ID string to play with
        runid = "";
        basecaller = "";
        start_time = "";
        if (tags.RG != NULL) {
            rg_info = create_rg_info(tags.RG);
            if (rg_info->runid != NULL) {
                runid = rg_info->runid;
            }
            if (rg_info->basecaller != NULL) {
                basecaller = rg_info->basecaller;
            }
        } else if (tags.RD != NULL) {
            runid = tags.RD;
        }

        if (tags.st != NULL) {
            start_time = tags.st;
        }
        kh_counter_increment(runids, runid);
        kh_counter_increment(basecallers, basecaller);

        // write a record for unmapped/unplaced
        if (b->core.flag & BAM_FUNMAP) {
            if (unmapped) {
                // an unmapped read can still have a RNAME and POS, but we
                // ignore that here, because its not a thing we care about
                char* qname = bam_get_qname(b);
                uint32_t read_length = b->core.l_qseq;
                float mean_quality = mean_qual_from_bam(bam_get_qual(b), read_length);
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
                        tags.dx
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
                        tags.dx
                    );
                }
                // add to flagstat counts if required
                if (flag_counts != NULL) {
                    process_flagstat_counts(b, flag_counts->unmapped, tags.dx);
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
            process_flagstat_counts(b, counts, tags.dx);
        }

        // only take "good" primary alignments for further processing
        if (b->core.flag & (NOTPRIMARY | BAM_FQCFAIL | BAM_FDUP)) {
            goto FINISH_READ;
        }
        char* qname = bam_get_qname(b);

        size_t* stats = create_cigar_stats(b);
        size_t match, ins, delt;
        // some aligners like to get fancy
        match = stats[BAM_CMATCH] + stats[BAM_CEQUAL] + stats[BAM_CDIFF];
        ins = stats[BAM_CINS];
        delt = stats[BAM_CDEL];
        size_t sub = tags.NM - ins - delt;
        size_t length = match + ins + delt;
        float iden = 100 * ((float)(match - sub)) / match;
        float acc = 100 - 100 * ((float)(tags.NM)) / length;
        // some things we've seen go wrong
        // explode now because there is almost certainly something wrong with the tags
        // and calling add_qual_count with a value less than zero will cause a segfault
        if (iden < 0.0 || acc < 0.0 || (size_t)tags.NM > length) {
            fprintf(stderr, "Read '%s' appears to contain implausible alignment information\n", qname);
            exit(EXIT_FAILURE);
        }
        // we only deal in primary/soft-clipped alignments so length
        // of qseq member is the length of the intact query sequence.
        uint32_t read_length = b->core.l_qseq;
        size_t qstart = get_query_start(b);
        size_t qend = get_query_end(b);
        // get mean quality score, from tag or recompute
        float mean_quality = tags.qs;
        if (mean_quality < 0 || force_recalc_qual) {
            mean_quality = mean_qual_from_bam_naive(bam_get_qual(b), read_length);
        }

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
                if (tags.pi == -1 && tags.pt >= 0) {
                    polya_len = tags.pt;
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
                match, ins, delt, sub, iden, acc, tags.dx);
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
                match, ins, delt, sub, iden, acc, tags.dx);
        }
		free(stats);

FINISH_READ:
        destroy_rg_info(rg_info);
        rg_info = NULL;
        runid = NULL;
        basecaller = NULL;
        start_time = NULL;
        free_bam_tags(&tags);
    }

    destroy_bam_iter_data(bam);
    bam_destroy1(b);

    return;
}
