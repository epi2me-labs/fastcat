#define _GNU_SOURCE
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/thread_pool.h"

#include "bamiter.h"
#include "common.h"
#include "counts.h"
#include "args.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table



/** Generates base counts from a region of a bam.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @returns void. Prints output to stdout.
 *
 */
void calculate_pileup(
        const char **bam_file, const char *chr, int start, int end,
        const char *read_group, const char tag_name[2], const int tag_value) {

    // setup bam reading - reuse our pileup structure, but actually just need iterator
    int nfile = 0; for (; bam_file[nfile]; nfile++);
    mplp_data **data = xalloc(nfile, sizeof(mplp_data*), "bam files");
    for (size_t i = 0; i < nfile; ++i) {
        data[i] = create_bam_iter_data(
            (const char *) bam_file[i], chr, start, end, read_group, tag_name, tag_value);
        if (data[i] == NULL) return NULL;
    }

    const bam_pileup1_t **plp = xalloc(nfile, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid;

    bam1_t *b = bam_init1();

    for (int fnumber=0; fnumber<nfile; fnumber++) {
        mplp_data bam = data[fnumber];
        while ((res = bam_itr_next(bam->fp, bam->iter, b) >= 0)) {
            const char *c_name = data[0]->hdr->target_name[tid];
            if (strcmp(c_name, chr) != 0) continue;
            if (pos < start) continue;
            if (pos >= end) break;
            fprintf(stdout, "%s\n", b->qname);
        }
    }

    for (size_t i = 0; i < nfile; ++i) {
        destroy_bam_iter_data(data[i]);
    }
    free(data);

    return;
}


typedef struct twarg {
    arguments_t args;
    const char *chr;
    int start;
    int end;
} twarg;


void *pileup_worker(void *arg) {
    twarg j = *(twarg *)arg;
    process_bams(
        j.args.bam, j.chr, j.start, j.end,
        j.args.read_group, j.args.tag_name, j.args.tag_value,
        j.args.lowthreshold, j.args.highthreshold, j.args.mod_base.code);
    free(arg);
    return pileup;
}


/* Process and print a single region using a threadpool
 *
 * @param args program arguments.
 * @param chr reference sequence to process.
 * @param start reference coordinate to process (0-based).
 * @param end reference coordiate to process (exclusive).
 * @param ref reference sequence.
 *
 */
#ifdef NOTHREADS
void process_region(arguments_t args, const char *chr, int start, int end, char *ref) {
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    process_bams(
        args.bam, chr, start, end,
        args.read_group, args.tag_name, args.tag_value,
        args.lowthreshold, args.highthreshold, args.mod_base.code);
}
#else
void process_region(arguments_t args, const char *chr, int start, int end, char *ref) {
    //TODO: Implement this.
    fprintf(stderr, "Called unimplemented threaded implementation.\n");
    exit(EXIT_FAILURE);
    //fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    //// create thread pool
    //hts_tpool *p = hts_tpool_init(args.threads);
    //hts_tpool_process *q = hts_tpool_process_init(p, 2 * args.threads, 0);
    //hts_tpool_result *r;
    //const int width = 1000000;

    //int nregs = 1 + (end - start) / width; float done = 0;
    //for (int rstart = start; rstart < end; rstart += width) {
    //    twarg *tw_args = xalloc(1, sizeof(*tw_args), "thread worker args");  // freed in worker
    //    tw_args->args = args;
    //    tw_args->chr = chr; tw_args->start = rstart; tw_args->end=min(rstart + width, end);
    //    int blk;
    //    do {
    //        blk = hts_tpool_dispatch2(p, q, pileup_worker, tw_args, 1);
    //        if ((r = hts_tpool_next_result(q))) {
    //            plp_data res = (plp_data)hts_tpool_result_data(r);
    //            if (res != NULL) {
    //                print_bedmethyl(
    //                    res, ref, 0,
    //                    args.extended, args.mod_base.abbrev, args.mod_base.base, args.cpg);
    //                destroy_plp_data(res);
    //                done++;
    //                fprintf(stderr, "\r%.1f %%", 100*done/nregs);
    //            }
    //            hts_tpool_delete_result(r, 0);
    //        }
    //    } while (blk == -1);
    //}

    //// wait for jobs, then collect.
    //hts_tpool_process_flush(q);
    //while ((r = hts_tpool_next_result(q))) {
    //    plp_data res = (plp_data)hts_tpool_result_data(r);
    //    if (res != NULL) {
    //        print_bedmethyl(
    //            res, ref, 0,
    //            args.extended, args.mod_base.abbrev, args.mod_base.base, args.cpg);
    //        destroy_plp_data(res);
    //        done++;
    //        fprintf(stderr, "\r%.1f %%", 100*done/nregs);
    //    }
    //    hts_tpool_delete_result(r, 0);
    //}
    //fprintf(stderr, "\r100 %%  ");
    //fprintf(stderr, "\n");
    //// clean up pool
    //hts_tpool_process_destroy(q);
    //hts_tpool_destroy(p);
}
#endif

