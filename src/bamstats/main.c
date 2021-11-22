// bamstats program

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"

#include "readstats.h"
#include "args.h"
#include "common.h"


void write_header() {
    fprintf(stdout,
"name\tref\tcoverage\tref_coverage\t"\
"\tqstart\tqend\trstart\trend\t"\
"aligned_ref_len\tdirection\tlength\tread_length"\
"match\tins\tdel\tsub\tiden\tacc\n");
}


int main(int argc, char *argv[]) {
    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
#ifdef NOTHREADS
    if (args.threads != 1) {
        fprintf(
            stderr,
            "--threads set to %d, but threading not supported by this build.\n", args.threads);
    }
#endif

    // large basecaller runs can produce more files than a single
    // process can open, check this ahead of time.
#ifndef WASM
    struct rlimit reslimit;
    int nfile = 0; for (; args.bam[nfile]; nfile++);
    if (getrlimit(RLIMIT_NOFILE, &reslimit) == 0) {
        if (nfile * args.threads > reslimit.rlim_cur - 100) {
            fprintf(stderr,
                "ERROR: Too many BAM files provided (%i). Try running "
                "samtools merge on subsets of files to produce fewer files", nfile);
            exit(EXIT_FAILURE);
        }
    }
#endif

    if (nfile > 1) {
        fprintf(stderr, "WARNING: Results from multiple files will not be coordinate sorted.\n");
    }

    write_header();

    htsFile *fps[args.threads];
    hts_idx_t *idxs[args.threads];
    sam_hdr_t *hdrs[args.threads];
    for (size_t i=0; i<args.threads; ++i) {
        fps[i] = hts_open(args.bam[0], "rb");
        idxs[i] = sam_index_load(fps[i], args.bam[0]);
        hdrs[i] = sam_hdr_read(fps[i]);
        if (hdrs[i] == 0 || idxs[i] == 0 || fps[i] == 0) {
            fprintf(stderr, "Failed to read .bam file '%s'.\n", args.bam[0]);
            exit(EXIT_FAILURE);
        }
    }


    if (args.region == NULL) {
        // process all regions
        for (int i=0; i < hdrs[0]->n_targets; ++i) {
            const char* chr = sam_hdr_tid2name(hdrs[0], i);
            size_t ref_length = (size_t)sam_hdr_tid2len(hdrs[0], i);
            process_region(
                fps, idxs, hdrs,
                args, chr, 0, ref_length, NULL);
        }
    } else {
        // process given region
        int start, end;
        char *chr = xalloc(strlen(args.region) + 1, sizeof(char), "chr");
        strcpy(chr, args.region);
        char *reg_chr = (char *) hts_parse_reg(chr, &start, &end);
        // start and end now zero-based end exclusive
        if (reg_chr) {
            *reg_chr = '\0';  // sets chr to be terminated at correct point
        } else {
            fprintf(stderr, "ERROR: Failed to parse region: '%s'.\n", args.region);
            exit(EXIT_FAILURE);
        }
        int tid = sam_hdr_name2tid(hdrs[0], chr);
        if (tid < 0) {
            fprintf(stderr, "ERROR: Failed to find reference '%s' in BAM header.\n", chr);
            exit(EXIT_FAILURE);
        }
        size_t ref_length = (size_t)sam_hdr_tid2len(hdrs[0], tid);
        end = min(end, ref_length);
        process_region(
            fps, idxs, hdrs,    
            args, chr, start, end, NULL);
        free(chr);
    }

    for (size_t i=0; i<args.threads; ++i) {
        hts_close(fps[i]);
        hts_idx_destroy(idxs[i]);
        sam_hdr_destroy(hdrs[i]);
    }

    clock_t end = clock();
    fprintf(stderr, "Total time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
