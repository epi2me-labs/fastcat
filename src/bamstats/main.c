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
#include "regiter.h"
#include "../fc_log.h"
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
    fc_log_config log_cfg = {
        .module = "BAMSTATS",
        .level = args.log_level,
        .target = args.log_file,
        .color = args.log_mono ? FC_LOG_COLOR_OFF : FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);
#ifdef NOTHREADS
    if (args.threads != 1) {
        fprintf(
            stderr,
            "--threads set to %d, but threading not supported by this build.\n", args.threads);
    }
#endif

    // set up output root (must not pre-exist)
    int can_output = can_make_dir(args.output_dir);
    if (can_output != 1) {
        FC_LOG_ERROR("Output directory '%s' already exists or cannot be created.", args.output_dir);
        exit(EXIT_FAILURE);
    }
    if (mkdir_hier(args.output_dir) == -1) {
        FC_LOG_ERROR("Failed to create output directory '%s'.", args.output_dir);
        exit(EXIT_FAILURE);
    }


    htsFile *fp = hts_open(args.bam, "rb");
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || fp == 0) {
        FC_LOG_ERROR("Failed to read .bam file '%s'.", args.bam);
        exit(EXIT_FAILURE);
    }

    htsThreadPool p = {NULL, 0};
    if (args.threads > 1) {
        FC_LOG_INFO("Using %d threads", args.threads);
        p.pool = hts_tpool_init(args.threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    fc_stats_engine_config stats_cfg = {
        .track_alignment = true,
        .track_polya = args.poly_a,
        .polya_cover = args.poly_a_cover,
        .polya_qual = args.poly_a_qual,
        .polya_rev = args.poly_a_rev
    };
    fc_writer *writer = fc_writer_create(args.demultiplex, &stats_cfg);
    if (writer == NULL) exit(EXIT_FAILURE);

    bool region_mode = (args.region != NULL) || (args.bed != NULL);
    if (args.per_read) {
        fc_perread_sink_opts perread_opts = {
            .output_root = args.output_dir,
            .sample = args.sample,
            .include_alignment = true,
            .use_barcodes = args.demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    }
    // per-file summary (single entry for bamstats)
    fc_writer_add_sink(writer, fc_perfile_sink_create(args.output_dir, args.sample));
    fc_hist_sink_opts hist_opts = {
        .output_root = args.output_dir,
        .use_barcodes = args.demultiplex,
        .unmapped_override = false,
        .write_alignment_metrics = true,
        .write_polya = args.poly_a,
        .engines = writer->stats_engines,
        .engine_count = writer->stats_engine_count
    };
    fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));
    fc_counter_sink_opts counter_opts = {
        .output_root = args.output_dir,
        .sample = args.sample,
        .use_barcodes = args.demultiplex,
        .barcode_count = writer->barcode_count
    };
    fc_writer_add_sink(writer, fc_counter_sink_create(&counter_opts));
    if (args.coverage) {
        fc_coverage_sink_opts cov_opts = {
            .output_root = args.output_dir,
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
            .n_segments = args.n_segments,
            .use_barcodes = args.demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_coverage_sink_create(&cov_opts));
    }
    fc_flagstats_sink_opts flag_opts = {
        .region_mode = region_mode,
        .output_root = args.output_dir,
        .sample = args.sample,
        .hdr = hdr,
        .use_barcodes = args.demultiplex,
        .barcode_count = writer->barcode_count
    };
    fc_writer_add_sink(writer, fc_flagstats_sink_create(&flag_opts));

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
        // proces given region / BED
        hts_idx_t *idx = sam_index_load(fp, args.bam);
        if (idx == 0){
            FC_LOG_ERROR("Cannot find index file for '%s', which is required for processing by region.", args.bam);
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
        FC_LOG_INFO("Processed %d regions", rit.n_regions);

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
    FC_LOG_INFO("Processed %" PRIu64 " records in 1 files.", total_records);
    FC_LOG_INFO("Parsing/filtering summary:");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        FC_LOG_INFO("%s\t%" PRIu64, failure_type[i], file_stats.failures[i]);
    }

    fc_writer_destroy(writer);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    destroy_args(&args);

    clock_t end = clock();
    FC_LOG_INFO("Total CPU time: %fs", (double)(end - begin) / CLOCKS_PER_SEC);
    fc_log_destroy();
    exit(EXIT_SUCCESS);
}
