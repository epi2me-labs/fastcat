// bamstats program

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "args.h"
#include "common.h"
#include "readstats.h"
#include "regiter.h"
#include "../fc_iter.h"
#include "../fc_sinks.h"
#include "../fc_writer.h"

static void init_bam_file_stats(fc_file_stats *stats) {
    if (!stats) return;
    memset(stats, 0, sizeof(*stats));
    stats->min_length = SIZE_MAX;
}

static void merge_bam_file_stats(fc_file_stats *dst, const fc_file_stats *src,
        double *qsum_total, size_t *n_total) {
    if (!dst || !src) return;
    dst->n_seqs += src->n_seqs;
    dst->n_bases += src->n_bases;
    if (src->n_seqs > 0) {
        if (dst->min_length == SIZE_MAX) {
            dst->min_length = src->min_length;
        } else {
            dst->min_length = min(dst->min_length, src->min_length);
        }
        dst->max_length = max(dst->max_length, src->max_length);
    }
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        dst->failures[i] += src->failures[i];
    }
    if (qsum_total && n_total && src->n_seqs > 0) {
        *qsum_total += src->mean_quality * (double)src->n_seqs;
        *n_total += src->n_seqs;
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

    htsFile *fp = hts_open(args.bam, "rb");
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || fp == 0) {
        fprintf(stderr, "ERROR: Failed to read .bam file '%s'.\n", args.bam);
        exit(EXIT_FAILURE);
    }

    htsThreadPool p = {NULL, 0};
    if (args.threads > 1) {
        fprintf(stderr, "Using %d threads\n", args.threads);
        p.pool = hts_tpool_init(args.threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    flag_stats* flag_counts = NULL;
    if (args.flagstats != NULL) {
        flag_counts = create_flag_stats(
            args.region == NULL ? hdr->n_targets : 1, args.unmapped
        );
    }

    fc_writer *writer = fc_writer_create();
    if (writer == NULL) exit(EXIT_FAILURE);
    fc_writer_set_barcode_count(writer, 1, false);

    fc_stats_engine_config stats_cfg = {
        .track_alignment = true,
        .track_unmapped = args.unmapped,
        .track_polya = args.poly_a,
        .polya_cover = args.poly_a_cover,
        .polya_qual = args.poly_a_qual,
        .polya_rev = args.poly_a_rev,
        .create_default = true
    };
    fc_writer_configure_stats(writer, &stats_cfg, false);

    bool region_mode = (args.region != NULL) || (args.bed != NULL);
    fc_perread_sink_opts perread_opts = {
        .path = NULL,
        .sample = args.sample,
        .include_alignment = true,
        .output_unmapped = args.unmapped
    };
    fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    if (args.histograms) {
        size_t engine_count = 0;
        fc_stats_engine **stats_engines = fc_writer_get_stats_engines(writer, &engine_count);
        fc_hist_sink_opts hist_opts = {
            .hist_dir = args.histograms,
            .demux_dir = NULL,
            .use_barcodes = false,
            .create_hist_dir = false,
            .write_unmapped = args.unmapped && !region_mode,
            .write_alignment_metrics = true,
            .write_polya = args.poly_a,
            .unmapped_is_primary = false,
            .engines = stats_engines,
            .engine_count = engine_count
        };
        fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));
    }
    if (args.runids || args.basecallers) {
        fc_writer_add_sink(writer, fc_counter_sink_create(args.runids, args.basecallers, args.sample));
    }
    if (args.coverage) {
        fc_coverage_sink_opts cov_opts = {
            .out_dir = args.coverages,
            .per_base = true,
            .use_cigar = false,
            .exclude_flags = -1,
            .include_flags = -1,
            .by_strand = true,
            .hdr = hdr,
            .pool = &p,
            .bed_files = args.coverage_beds,
            .bed_names = args.coverage_names,
            .n_beds = args.n_coverage_beds,
            .thresholds = args.coverage_thresholds,
            .n_thresholds = args.n_coverage_thresholds,
            .segments = args.segments,
            .n_segments = args.n_segments
        };
        fc_writer_add_sink(writer, fc_coverage_sink_create(&cov_opts));
    }
    if (flag_counts != NULL) {
        fc_flagstats_sink_opts flag_opts = {
            .stats = flag_counts,
            .region_mode = region_mode,
            .track_unmapped = args.unmapped,
            .out_path = args.flagstats,
            .sample = args.sample,
            .hdr = hdr
        };
        fc_writer_add_sink(writer, fc_flagstats_sink_create(&flag_opts));
    }

    fc_bam_iter_opts iter_opts = {
        .read_group = args.read_group,
        .tag_name = { args.tag_name[0], args.tag_name[1] },
        .tag_value = args.tag_value,
        .force_recalc_qual = args.force_recalc_qual,
        .min_length = 0,
        .max_length = (size_t)-1,
        .min_qscore = 0.0f,
        .dust = false,
        .dust_t = 0,
        .dust_w = 0,
        .max_dust = 0.0f
    };

    fc_file_stats file_stats = {0};
    init_bam_file_stats(&file_stats);
    double qsum_total = 0.0;
    size_t n_total = 0;
    bool saw_stream_error = false;
    bool saw_unknown_error = false;

    fc_writer_begin_file(writer, args.bam);

    if (!region_mode) {
        fc_bam_iter iter = {0};
        if (fc_bam_iter_open(&iter, fp, NULL, hdr, NULL, 0, INT64_MAX, &iter_opts) == 0) {
            fc_read *fr = NULL;
            while ((fr = fc_bam_iter_next(&iter, hdr)) != NULL) {
                fc_writer_add_read(writer, fr, args.bam, hdr);
                fc_read_destroy(fr);
            }
            fc_file_stats iter_stats = {0};
            fc_bam_iter_finish(&iter, &iter_stats);
            saw_stream_error |= (iter_stats.failures[F_STREAM_ERROR] > 0);
            saw_unknown_error |= (iter_stats.failures[F_UNKNOWN_ERROR] > 0);
            iter_stats.failures[F_FILE_OK] = 0;
            iter_stats.failures[F_STREAM_ERROR] = 0;
            iter_stats.failures[F_UNKNOWN_ERROR] = 0;
            merge_bam_file_stats(&file_stats, &iter_stats, &qsum_total, &n_total);
            fc_bam_iter_close(&iter);
        } else {
            saw_unknown_error = true;
        }

    } else {
        hts_idx_t *idx = sam_index_load(fp, args.bam);
        if (idx == 0){
            fprintf(stderr, "ERROR: Cannot find index file for '%s', which is required for processing by region.\n", args.bam);
            exit(EXIT_FAILURE);
        }

        regiter rit = init_region_iterator(args.bed, args.region, hdr);
        int check = 0;
        while ((check = next_region(&rit)) != -1) {
            if (check == -2 && args.bed == NULL) {
                exit(EXIT_FAILURE);
            }
            if (check != 0) continue;

            fc_bam_iter iter = {0};
            if (fc_bam_iter_open(&iter, fp, idx, hdr, rit.chr, rit.start, rit.end, &iter_opts) == 0) {
                fc_read *fr = NULL;
                while ((fr = fc_bam_iter_next(&iter, hdr)) != NULL) {
                    fc_writer_add_read(writer, fr, args.bam, hdr);
                    fc_read_destroy(fr);
                }
                fc_file_stats iter_stats = {0};
                fc_bam_iter_finish(&iter, &iter_stats);
                saw_stream_error |= (iter_stats.failures[F_STREAM_ERROR] > 0);
                saw_unknown_error |= (iter_stats.failures[F_UNKNOWN_ERROR] > 0);
                iter_stats.failures[F_FILE_OK] = 0;
                iter_stats.failures[F_STREAM_ERROR] = 0;
                iter_stats.failures[F_UNKNOWN_ERROR] = 0;
                merge_bam_file_stats(&file_stats, &iter_stats, &qsum_total, &n_total);
                fc_bam_iter_close(&iter);
            } else {
                saw_unknown_error = true;
            }
        }
        fprintf(stderr, "Processed %d regions\n", rit.n_regions);

        destroy_region_iterator(&rit);
        hts_idx_destroy(idx);
    }

    if (n_total > 0) {
        file_stats.mean_quality = qsum_total / (double)n_total;
    } else {
        file_stats.mean_quality = 0.0;
    }
    if (file_stats.min_length == SIZE_MAX) {
        file_stats.min_length = 0;
    }
    file_stats.failures[F_STREAM_ERROR] = saw_stream_error ? 1 : 0;
    file_stats.failures[F_UNKNOWN_ERROR] = saw_unknown_error ? 1 : 0;
    file_stats.failures[F_FILE_OK] = (!saw_stream_error && !saw_unknown_error) ? 1 : 0;

    fc_writer_end_file(writer, args.bam, &file_stats);
    fc_writer_finish(writer);

    uint64_t total_records =
        file_stats.failures[R_RECORD_OK]
        + file_stats.failures[R_TOO_LONG]
        + file_stats.failures[R_TOO_SHORT]
        + file_stats.failures[R_LOW_QUALITY]
        + file_stats.failures[R_DUST_MASKED];
    fprintf(stderr, "INFO   : Processed %" PRIu64 " records in 1 files.\n", total_records);
    fprintf(stderr, "\nParsing/filtering summary:\n");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        fprintf(stderr, "%s\t%" PRIu64 "\n", failure_type[i], file_stats.failures[i]);
    }

    fc_writer_destroy(writer);

    if (flag_counts != NULL) destroy_flag_stats(flag_counts);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    destroy_args(&args);

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
