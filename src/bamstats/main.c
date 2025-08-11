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

#include "args.h"
#include "common.h"
#include "readstats.h"
#include "regiter.h"


void write_header(const char* sample) {
    char *sn = sample == NULL ? "" : "\tsample_name";
    fprintf(stdout, 
        "name\trunid%s\tref\tcoverage\tref_coverage\t"\
        "qstart\tqend\trstart\trend\t"\
        "aligned_ref_len\tdirection\tlength\tread_length\tmean_quality\tstart_time\t"\
        "match\tins\tdel\tsub\tiden\tacc\tduplex\n",
        sn);
}

// stats array should have 8 entries
// total, primary, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP, BAM_FQCFAIL, BAM_FDUP, unused
// note: HTS spec makes a distinction between "unmapped" (flag & 4) and "unplaced". Unplaced
//       are not necessarily unmapped but lack definitive coords, this is mainly for paired-end
//       but we'll keep the distinction here.
void write_stats_header(FILE* fh, const char* sample) {
    char *sn = sample == NULL ? "" : "\tsample_name";
    fprintf(fh, "ref%s\ttotal\tprimary\tsecondary\tsupplementary\tunmapped\tqcfail\tduplicate\tduplex\tduplex_forming\n", sn);
}

static inline void write_stats(size_t *stats, const char* chr, const char* sample, FILE* fh) {
    if (fh != NULL) {
        if (sample == NULL) {
            fprintf(fh,
                "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
                chr, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7], stats[8]
            );
        } else {
            fprintf(fh,
                "%s\t%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
                chr, sample, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7], stats[8]
            );
        }
    }
}

static inline void write_counter(const char* fname, kh_counter_t *counter, const char* sample, const char* bam_fname, const char* column_name) {
    ensure_parent_dir_exists(fname);
    FILE* stats_fp = fopen(fname, "w");
    if (stats_fp == NULL) {
        fprintf(stderr, "ERROR: Cannot open file '%s' for writing.\n", fname);
        exit(EXIT_FAILURE);
    }
    fprintf(stats_fp, "filename\t");
    if (sample != NULL) fprintf(stats_fp, "sample_name\t");
    fprintf(stats_fp, "%s\tcount\n", column_name);
    for (khiter_t k = 0; k < kh_end(counter); ++k) {
        if (kh_exist(counter, k)) {
            fprintf(stats_fp, "%s\t", bam_fname);
            if (sample != NULL) fprintf(stats_fp, "%s\t", sample);
            fprintf(stats_fp, "%s\t%d\n", kh_key(counter, k), kh_val(counter, k));
        }
    }
    fclose(stats_fp);
}


void write_hist_stats(read_stats* stats, char* prefix, char* name) {
    char* path = calloc(strlen(prefix) + strlen(name) + 2, sizeof(char));
    sprintf(path, "%s/%s", prefix, name);
    ensure_parent_dir_exists(path);
    FILE* fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Cannot open file '%s' for writing.\n", path);
        exit(EXIT_FAILURE);
    }
    print_stats(stats, false, true, fp);
    fclose(fp); free(path);
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
    size_t nfile = 0; for (; args.bam[nfile]; nfile++);
    if (getrlimit(RLIMIT_NOFILE, &reslimit) == 0) {
        if (nfile * args.threads > reslimit.rlim_cur - 100) {
            fprintf(stderr,
                "ERROR: Too many BAM files provided (%zu). Try running "
                "samtools merge on subsets of files to produce fewer files", nfile);
            exit(EXIT_FAILURE);
        }
    }
#endif

    // Check all outputs do not pre-exist to avoid overwriting 
    int rtn = mkdir_hier(args.histograms);
    if (rtn == -1) {
        fprintf(stderr,
           "ERROR: Cannot create output directory '%s'. Check location is writeable and directory does not exist.\n",
           args.histograms);
        exit(EXIT_FAILURE);
    }

    if (args.flagstats != NULL && file_exists(args.flagstats)) {
        fprintf(stderr, "ERROR: Output file '%s' already exists, please remove it or use a different name.\n", args.flagstats);
        exit(EXIT_FAILURE);
    }
    if (args.runids != NULL && file_exists(args.runids)) {
        fprintf(stderr, "ERROR: Output file '%s' already exists, please remove it or use a different name.\n", args.runids);
        exit(EXIT_FAILURE);
    }
    if (args.basecallers != NULL && file_exists(args.basecallers)) {
        fprintf(stderr, "ERROR: Output file '%s' already exists, please remove it or use a different name.\n", args.basecallers);
        exit(EXIT_FAILURE);
    }

    // TODO: don't be lazy
    if (nfile > 1) {
        fprintf(stderr, "ERROR: Multiple input files detected, this program currently supports only a single file.\n");
        exit(EXIT_FAILURE);
    }

    write_header(args.sample);

    htsFile *fp = hts_open(args.bam[0], "rb");
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || fp == 0) {
        fprintf(stderr, "ERROR: Failed to read .bam file '%s'.\n", args.bam[0]);
        exit(EXIT_FAILURE);
    }

    htsThreadPool p = {NULL, 0};
    if (args.threads > 1 ) {
        fprintf(stderr, "Using %d threads\n", args.threads);
        p.pool = hts_tpool_init(args.threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    FILE* flagstats = NULL;
    flag_stats* flag_counts = NULL;
    if (args.flagstats != NULL) {
        ensure_parent_dir_exists(args.flagstats);
        flagstats = fopen(args.flagstats, "w");
        if (flagstats == NULL) {
            fprintf(stderr, "ERROR: Cannot open file '%s' for writing.\n", args.flagstats);
            exit(EXIT_FAILURE);
        }
        write_stats_header(flagstats, args.sample);
        flag_counts = create_flag_stats(
            args.region == NULL ? hdr->n_targets : 1, args.unmapped
        );
    }

    kh_counter_t *run_ids = kh_counter_init();
    kh_counter_t *basecallers = kh_counter_init();
    read_stats* length_stats = create_length_stats();
    read_stats* qual_stats = create_qual_stats(QUAL_HIST_WIDTH);
    read_stats* acc_stats = create_qual_stats(ACC_HIST_WIDTH);
    read_stats* cov_stats = create_qual_stats(COV_HIST_WIDTH);
    read_stats* polya_stats = args.poly_a ? create_length_stats() : NULL;

    // Prepare also for the unmapped reads
    read_stats* length_stats_unmapped = create_length_stats();
    read_stats* qual_stats_unmapped = create_qual_stats(QUAL_HIST_WIDTH);

    if (args.region == NULL && args.bed == NULL) {
        // iterate over the entire file
        process_bams(
            fp, NULL, hdr, args.sample,
            NULL, 0, INT64_MAX, true,
            args.read_group, args.tag_name, args.tag_value,
            flag_counts, args.unmapped,
            length_stats, qual_stats, acc_stats, cov_stats,
            length_stats_unmapped, qual_stats_unmapped,
            polya_stats, args.poly_a_cover, args.poly_a_qual, args.poly_a_rev,
            run_ids, basecallers, args.force_recalc_qual);

        // write flagstat counts if requested
        if (flag_counts != NULL) {
            for (int i=0; i < hdr->n_targets; ++i) {
                const char* chr = sam_hdr_tid2name(hdr, i);
                write_stats(flag_counts->counts[i], chr, args.sample, flagstats);
            }
            if (args.unmapped) {
                write_stats(flag_counts->unmapped, "*", args.sample, flagstats);
            }
        }
    } else {
        // process given region / BED
        hts_idx_t *idx = sam_index_load(fp, args.bam[0]);
        if (idx == 0){
            fprintf(stderr, "ERROR: Cannot find index file for '%s', which is required for processing by region.\n", args.bam[0]);
            exit(EXIT_FAILURE);
        }

        regiter rit = init_region_iterator(args.bed, args.region, hdr);
        int check = 0;
        while ((check = next_region(&rit)) != -1) {
            if (check == -2 && args.bed == NULL) {
                // we were given only a region, not a bed, and that region was garbage
                // => user error, should stop immediately
                exit(EXIT_FAILURE);
            }
            if (check != 0) continue;  // skip other errors

            process_bams(
                fp, idx, hdr, args.sample,
                rit.chr, rit.start, rit.end, true,
                args.read_group, args.tag_name, args.tag_value,
                flag_counts, args.unmapped,
                length_stats, qual_stats, acc_stats, cov_stats,
                length_stats_unmapped, qual_stats_unmapped,
                polya_stats, args.poly_a_cover, args.poly_a_qual, args.poly_a_rev,
                run_ids, basecallers, args.force_recalc_qual);
            if (flag_counts != NULL) {
                // TODO: regions might not be whole chromosomes...
                write_stats(flag_counts->counts[0], rit.chr, args.sample, flagstats);
            }
        }
        fprintf(stderr, "Processed %d regions\n", rit.n_regions);
        
        destroy_region_iterator(&rit);
        hts_idx_destroy(idx);
    }

    write_hist_stats(length_stats, args.histograms, "length.hist");
    write_hist_stats(qual_stats, args.histograms, "quality.hist");
    write_hist_stats(acc_stats, args.histograms, "accuracy.hist");
    write_hist_stats(cov_stats, args.histograms, "coverage.hist");
    if (polya_stats != NULL) {
        write_hist_stats(polya_stats, args.histograms, "polya.hist");
    } 

    // Save also histograms for the unmapped reads if requested
    // and if the user is not asking for a region
    if (args.unmapped && args.region == NULL){
        write_hist_stats(length_stats_unmapped, args.histograms, "length.unmap.hist");
        write_hist_stats(qual_stats_unmapped, args.histograms, "quality.unmap.hist");
    }

    // write runids summary
    if (args.runids != NULL) {
        write_counter(args.runids, run_ids, args.sample, args.bam[0], "run_id");
    } 
    // write basecallers summary
    if (args.basecallers != NULL) {
        write_counter(args.basecallers, basecallers, args.sample, args.bam[0], "basecaller");
    } 

    destroy_length_stats(length_stats);
    destroy_qual_stats(qual_stats);
    destroy_qual_stats(acc_stats);
    destroy_qual_stats(cov_stats);
    destroy_length_stats(length_stats_unmapped);
    destroy_qual_stats(qual_stats_unmapped);
    destroy_length_stats(polya_stats);
    kh_counter_destroy(basecallers);
    kh_counter_destroy(run_ids);

    if (flagstats != NULL) {
        fclose(flagstats);
    }

    if (flag_counts != NULL) destroy_flag_stats(flag_counts);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
