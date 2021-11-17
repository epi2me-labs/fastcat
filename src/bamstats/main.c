// bamstats program

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/faidx.h"

#include "readstats.h"
#include "args.h"
#include "common.h"


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

    // load ref sequence
    faidx_t *fai = fai_load(args.ref);
    if (fai == NULL) {
        fprintf(stderr,
            "ERROR: Failed to parse reference file\n");
        exit(EXIT_FAILURE);
    }
    if (args.region == NULL) {
        // process all regions
        int nseq = faidx_nseq(fai);
        for (int i = 0; i < nseq; ++i) {
            const char *chr = faidx_iseq(fai, i);
            int len = faidx_seq_len(fai, chr);
            int alen;
            char *ref = faidx_fetch_seq(fai, chr, 0, len, &alen);
            process_region(args, chr, 0, len, ref);
            free(ref);
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
        // simplify things for later by fetching whole chr
        int len;
        char *ref = fai_fetch(fai, chr, &len);
        if (len < 0) {
            fprintf(stderr, "ERROR: Failed to fetch reference region: '%s'.\n", args.region);
            exit(EXIT_FAILURE);
        }
        end = min(end, len);
        process_region(args, chr, start, end, ref);

        free(chr);
        free(ref);
    }
    fai_destroy(fai);
    clock_t end = clock();
    fprintf(stderr, "Total time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
