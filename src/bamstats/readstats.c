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
#include "htslib/thread_pool.h"

#include "../common.h"
#include "bamiter.h"
#include "readstats.h"
#include "args.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table

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


/** Generates alignment stats from a region of a bam.
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @returns void. Prints output to stdout.
 *
 */
void process_bams(
        htsFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
        const char *chr, int start, int end, bool overlap_start,
        const char *read_group, const char tag_name[2], const int tag_value) {

    // setup bam reading - reuse our pileup structure, but actually just need iterator
    mplp_data* bam = create_bam_iter_data(
        fp, idx, hdr,
        chr, start, end, overlap_start,
        read_group, tag_name, tag_value);
    if (bam == NULL) return;

    int res;
    bam1_t *b = bam_init1();
    uint8_t *tag;
    // find the target length for query below
    size_t ref_length = 0;
    int tid = sam_hdr_name2tid(hdr, chr);
    if (tid < 0) {
        fprintf(stderr, "Failed to find reference sequence '%s' in bam.\n", chr);
        exit(EXIT_FAILURE);
    }
    ref_length = (size_t)sam_hdr_tid2len(hdr, tid);

    while ((res = read_bam(bam, b) >= 0)) {
        char* qname = bam_get_qname(b);

        // get NM tag
        tag = bam_aux_get((const bam1_t*) b, "NM");
        if (tag == NULL){ // tag isn't present or is currupt
            fprintf(stderr, "Read '%s' does not contain 'NM' tag.\n", qname);
            exit(EXIT_FAILURE);
        }
        int NM = bam_aux2i(tag);
        if (errno == EINVAL) {
            fprintf(stderr, "Read '%s' contains non-integer 'NM' tag.\n", qname);
            exit(EXIT_FAILURE); 
        }
		size_t* stats = create_cigar_stats(b);
		size_t match, ins, delt;
        match = stats[0]; ins = stats[1]; delt = stats[2];
		size_t sub = NM - ins - delt;
		size_t length = match + ins + delt;
		float iden = 100 * ((float)(match - sub)) / match;
		float acc = 100 - 100 * ((float)(NM)) / length;
        // we only deal in primary/soft-clipped alignments so length
        // ok qseq member is the length of the intact query sequence.
        uint32_t read_length = b->core.l_qseq;
        size_t qstart = get_query_start(b);
        size_t qend = get_query_end(b);

        float coverage = 100 * ((float)(qend - qstart)) / read_length;
        size_t rstart = b->core.pos;
        size_t rend = bam_endpos(b);
        size_t aligned_ref_len = rend - rstart;
        float ref_cover = 100 * ((float)(aligned_ref_len)) / ref_length;
        char direction = "+-"[bam_is_rev(b)];

        fprintf(stdout,
            "%s\t%s\t%.4f\t%.4f\t" \
            "%lu\t%lu\t%lu\t%lu\t" \
            "%lu\t%c\t%lu\t%u\t" \
            "%lu\t%lu\t%lu\t%lu\t%.3f\t%.3f\n",
            qname, chr, coverage, ref_cover,
            qstart, qend, rstart, rend,
            aligned_ref_len, direction, length, read_length,
            match, ins, delt, sub, iden, acc);
		free(stats);
    }

    destroy_bam_iter_data(bam);
    bam_destroy1(b);

    return;
}

static int worker_id(hts_tpool *p) {
    int i;
    pthread_t s = pthread_self();
    for (i = 0; i < p->tsize; i++) {
        if (pthread_equal(s, p->t[i].tid))
            return i;
    }
    return -1;
}


typedef struct twarg {
    arguments_t args;
    const char *chr;
    int start;
    int end;
    bool overlap_start;
    hts_tpool *pool;
    htsFile **fps;
    hts_idx_t **idxs;
    sam_hdr_t **hdrs;
} twarg;


void *bam_worker(void *arg) {
    int *rtn = xalloc(1, sizeof(int), "worker_result");
    twarg j = *(twarg *)arg;
    int tid = worker_id(j.pool);
    if (tid < 0) {
        fprintf(stderr, "Could not retrieve thread ID in worker function\n");
        rtn[0] = 1; return rtn;
    }
    process_bams(
        j.fps[tid], j.idxs[tid], j.hdrs[tid],
        j.chr, j.start, j.end, j.overlap_start,
        j.args.read_group, j.args.tag_name, j.args.tag_value);
    free(arg);
    rtn[0] = 0;
    return rtn;
}


/* Process and print a single region using a threadpool
 *
 *  @param fp htsFile pointer
 *  @param idx hts_idx_t pointer
 *  @param hdr sam_hdr_t pointer
 *  @param args program arguments.
 *  @param chr reference sequence to process.
 *  @param start reference coordinate to process (0-based).
 *  @param end reference coordiate to process (exclusive).
 *  @param ref reference sequence.
 *
 */
#ifdef NOTHREADS
int process_region(
        htsFile **fps, hts_idx_t **idxs, sam_hdr_t **hdrs,
        arguments_t args, const char *chr, int start, int end, char *ref) {
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    process_bams(
        fps[0], idxs[0], hdrs[0],    
        chr, start, end,
        args.read_group, args.tag_name, args.tag_value);
    return 0;
}
#else
void process_region(
        htsFile **fps, hts_idx_t **idxs, sam_hdr_t **hdrs,
        arguments_t args, const char *chr, int start, int end, char *ref) {
    //TODO: Implement this properly
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    hts_tpool *p = hts_tpool_init(args.threads);
    hts_tpool_process *q = hts_tpool_process_init(p, 2 * args.threads, 0);
    hts_tpool_result *r;
    const int width = 100000;

    int nregs = 1 + (end - start) / width; float done = 0;
    bool overlap_start = true; // process reads overlapping start?
    for (int rstart = start; rstart < end; rstart += width) {
        twarg *tw_args = xalloc(1, sizeof(*tw_args), "thread worker args");  // freed in worker
        tw_args->args = args;
        tw_args->chr = chr; tw_args->start = rstart; tw_args->end=min(rstart + width, end);
        tw_args->overlap_start = overlap_start;
        tw_args->fps = fps; tw_args->hdrs = hdrs; tw_args->idxs = idxs;
        tw_args->pool = p;
        int blk;
        do {
            blk = hts_tpool_dispatch2(p, q, bam_worker, tw_args, 1);
            if ((r = hts_tpool_next_result(q))) {
                int* res = (int*)hts_tpool_result_data(r);
                if (*res == 0) {
                    done++;
                    fprintf(stderr, "\r%.1f %%", 100*done/nregs);
                } else {
                    fprintf(stderr, "\nPool worker failed with: %d\n", *res);
                }
                free(res);
                hts_tpool_delete_result(r, 0);
            }
        } while (blk == -1);
        overlap_start = false; // all subsequent blocks, don't double count reads
    }

    // wait for jobs, then collect.
    hts_tpool_process_flush(q);
    while ((r = hts_tpool_next_result(q))) {
        int* res = (int*)hts_tpool_result_data(r);
        if (*res == 0) {
            done++;
            fprintf(stderr, "\r%.1f %%", 100*done/nregs);
        } else {
            fprintf(stderr, "\nPool worker failed with: %d\n", *res);
        }
        free(res);
        hts_tpool_delete_result(r, 0);
    }
    fprintf(stderr, "\r100 %%  ");
    fprintf(stderr, "\n");
    // clean up pool
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
}
#endif

