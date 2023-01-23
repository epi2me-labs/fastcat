// bamstats program

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "readstats.h"
#include "args.h"
#include "common.h"


void write_header(const char* sample) {
    if (sample == NULL) {
        fprintf(stdout,
            "name\tref\tcoverage\tref_coverage\t"\
            "qstart\tqend\trstart\trend\t"\
            "aligned_ref_len\tdirection\tlength\tread_length\tmean_quality\t"\
            "match\tins\tdel\tsub\tiden\tacc\n");
    } else {
        fprintf(stdout,
            "name\tsample_name\tref\tcoverage\tref_coverage\t"\
            "qstart\tqend\trstart\trend\t"\
            "aligned_ref_len\tdirection\tlength\tread_length\tmean_quality\t"\
            "match\tins\tdel\tsub\tiden\tacc\n");
    }
}

// stats array should have 8 entries
// total, primary, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP, BAM_FQCFAIL, BAM_FDUP, unused
// note: HTS spec makes a distinction between "unmapped" (flag & 4) and "unplaced". Unplaced
//       are not necessarily unmapped but lack definitive coords, this is mainly for paired-end
//       but we'll keep the distinction here.
void write_stats_header(FILE* fh, const char* sample) {
    if (sample == NULL) {
        fprintf(fh, "ref\ttotal\tprimary\tsecondary\tsupplementary\tunmapped\tqcfail\tduplicate\n");
    } else {
        fprintf(fh, "ref\tsample_name\ttotal\tprimary\tsecondary\tsupplementary\tunmapped\tqcfail\tduplicate\n");
    }
}

static inline void write_stats(size_t *stats, const char* chr, const char* sample, FILE* fh) {
    if (fh != NULL) {
        if (sample == NULL) {
            fprintf(fh,
                "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
                chr, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6]
            );
        } else {
            fprintf(fh,
                "%s\t%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
                chr, sample, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6]
            );
        }
    }
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

    write_header(args.sample);

    htsFile *fp = hts_open(args.bam[0], "rb");
    hts_idx_t *idx = sam_index_load(fp, args.bam[0]);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || idx == 0 || fp == 0) {
        fprintf(stderr, "Failed to read .bam file '%s'.\n", args.bam[0]);
        exit(EXIT_FAILURE);
    }

    htsThreadPool p = {NULL, 0};
    if (args.threads > 1 ) {
        p.pool = hts_tpool_init(args.threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    FILE* flagstats = NULL;
    if (args.flagstats != NULL) {
        flagstats = fopen(args.flagstats, "w");
        write_stats_header(flagstats, args.sample);
    }

    size_t* flag_counts = xalloc(8, sizeof(size_t), "counts");
    if (args.region == NULL) {
        // process all regions
        for (int i=0; i < hdr->n_targets; ++i) {
            const char* chr = sam_hdr_tid2name(hdr, i);
            size_t ref_length = (size_t)sam_hdr_tid2len(hdr, i);
            process_bams(
                fp, idx, hdr, args.sample,
                chr, 0, ref_length, true,
                args.read_group, args.tag_name, args.tag_value,
                flag_counts, args.unmapped);
            write_stats(flag_counts, chr, args.sample, flagstats);
            memset(flag_counts, 0, 8 * sizeof(size_t));
        }
        // Also do unplaced reads
        if (args.unmapped) {
            process_bams(
                fp, idx, hdr, args.sample,
                "*", 0, INT64_MAX, true,
                args.read_group, args.tag_name, args.tag_value,
                flag_counts, args.unmapped);
            write_stats(flag_counts, "*", args.sample, flagstats);
            memset(flag_counts, 0, 8 * sizeof(size_t));
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
        int tid = sam_hdr_name2tid(hdr, chr);
        if (tid < 0) {
            fprintf(stderr, "ERROR: Failed to find reference '%s' in BAM header.\n", chr);
            exit(EXIT_FAILURE);
        }
        size_t ref_length = (size_t)sam_hdr_tid2len(hdr, tid);
        end = min(end, ref_length);
        process_bams(
            fp, idx, hdr, args.sample,
            chr, start, end, true,
            args.read_group, args.tag_name, args.tag_value,
            flag_counts, args.unmapped);
        write_stats(flag_counts, chr, args.sample, flagstats);
        memset(flag_counts, 0, 8 * sizeof(size_t));
        free(chr);
    }
    free(flag_counts);
    
    if (flagstats != NULL) {
        fclose(flagstats);
    }

    sam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
