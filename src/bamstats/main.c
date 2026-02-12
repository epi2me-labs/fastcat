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

    fc_iter_opts iter_opts = {
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
        .max_dust = 0.0f,
        .region = region_mode ? args.region : NULL,
        .bed = region_mode ? args.bed : NULL
    };

    // collect inputs via unified iterator, but enforce single BAM
    char *paths[2] = { (char *)args.bam, NULL };
    fc_iter file_iter = {0};
    if (fc_iter_init(&file_iter, paths, 1, 0, false) != 0 || file_iter.n_files == 0) {
        FC_LOG_ERROR("No BAM inputs found.");
        exit(EXIT_FAILURE);
    }
    if (file_iter.n_files != 1) {
        FC_LOG_ERROR("bamstats supports a single BAM input; found %zu inputs.", file_iter.n_files);
        exit(EXIT_FAILURE);
    }
    if (file_iter.files[0].kind != FC_FILE_BAM) {
        FC_LOG_ERROR("Input '%s' is not a BAM/SAM file.", file_iter.files[0].path);
        exit(EXIT_FAILURE);
    }

    fc_file_stats global_stats = {0};
    fc_file_stats_init(&global_stats);
    int status = EXIT_SUCCESS;

    const fc_file_entry *f = NULL;
    while (fc_iter_next_file(&file_iter, &f)) {
        fc_writer_begin_file(writer, f->path);

        if (fc_iter_open_record_cursor(&file_iter, &iter_opts) != 0) {
            status = EXIT_FAILURE;
            fc_iter_merge_file_stats(&file_iter, &global_stats);
            fc_writer_end_file(writer, f->path, fc_iter_stats(&file_iter));
            continue;
        }

        fc_read *fr = NULL;
        while ((fr = fc_iter_next_record(&file_iter, NULL)) != NULL) {
            fc_writer_add_read(writer, fr, f->path, file_iter.bam_hdr);
            fc_read_destroy(fr);
        }

        int fin_rc = fc_iter_finish_file(&file_iter);
        if (fin_rc != 0) status = EXIT_FAILURE;
        fc_iter_merge_file_stats(&file_iter, &global_stats);
        fc_writer_end_file(writer, f->path, fc_iter_stats(&file_iter));
    }

    fc_iter_close(&file_iter);

    fc_writer_finish(writer);

    uint64_t total_records =
        global_stats.failures[R_RECORD_OK]
        + global_stats.failures[R_TOO_LONG]
        + global_stats.failures[R_TOO_SHORT]
        + global_stats.failures[R_LOW_QUALITY]
        + global_stats.failures[R_DUST_MASKED];
    FC_LOG_INFO("Processed %" PRIu64 " records in 1 files.", total_records);
    if (status != EXIT_SUCCESS) {
        FC_LOG_WARN("Error processing files.");
    }
    FC_LOG_INFO("Parsing/filtering summary:");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        FC_LOG_INFO("%s\t%" PRIu64, failure_type[i], global_stats.failures[i]);
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
    exit(status);
}
