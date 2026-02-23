#include "fc_program.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"
#include "fc_iter.h"
#include "fc_log.h"
#include "fc_sinks.h"
#include "fc_writer.h"
#include "fc_filter.h"
#include "htslib/thread_pool.h"

static int process_current_file(
    fc_iter *file_iter,
    fc_writer *writer,
    const char *path,
    const sam_hdr_t *hdr_for_records,
    const fc_iter_opts *iter_opts,
    bool cursor_already_open,
    bool enforce_unmapped,
    bool *hard_error,
    fc_file_stats *global_stats) {
    int status = EXIT_SUCCESS;

    fc_writer_begin_file(writer, path);

    if (!cursor_already_open && fc_iter_open_record_cursor(file_iter, iter_opts) != 0) {
        status = EXIT_FAILURE;
        fc_iter_merge_file_stats(file_iter, global_stats);
        fc_writer_end_file(writer, path, fc_iter_stats(file_iter));
        return status;
    }

    fc_read *fr = NULL;
    while ((fr = fc_iter_next_record(file_iter)) != NULL) {
        if (enforce_unmapped && !fc_is_unmapped(fr->bam)) {
            FC_LOG_ERROR(
                "Input BAM contains mapped reads (e.g. '%s'). "
                "fastcat fastq with --input-fmt bam only accepts unaligned BAM input.",
                bam_get_qname(fr->bam));
            if (hard_error) *hard_error = true;
            status = EXIT_FAILURE;
            fc_read_destroy(fr);
            break;
        }
        fc_writer_add_read(writer, fr, path, hdr_for_records);
        fc_read_destroy(fr);
    }

    if (fc_iter_finish_file(file_iter) != 0) {
        status = EXIT_FAILURE;
    }
    fc_iter_merge_file_stats(file_iter, global_stats);
    fc_writer_end_file(writer, path, fc_iter_stats(file_iter));
    return status;
}

int fc_run_program(fc_program_kind program, const fc_catstats_args *args) {
    clock_t begin = clock();
    int status = EXIT_SUCCESS;

    fc_log_config log_cfg = {
        .module = (program == FC_PROGRAM_BAMSTATS) ? "BAMSTATS" : "FASTCAT",
        .level = args->log_level,
        .target = args->log_file,
        .color = args->log_mono ? FC_LOG_COLOR_OFF : FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    // initialised below if we successfully create outputs
    htsThreadPool pool = {NULL, 0};
    fc_writer* writer = NULL;
    fc_iter file_iter = {0};
    bool hard_error = false;

#ifdef NOTHREADS
    if (args->threads != 1) {
        fprintf(stderr, "--threads set to %d, but threading not supported by this build.\n", args->threads);
    }
#endif

    if (can_make_dir(args->output_dir) != 1) {
        FC_LOG_ERROR("Output directory '%s' already exists or cannot be created.", args->output_dir);
        status = EXIT_FAILURE;
        goto cleanup;
    }
    if (mkdir_hier(args->output_dir) == -1) {
        FC_LOG_ERROR("Failed to create output directory '%s'.", args->output_dir);
        status = EXIT_FAILURE;
        goto cleanup;
    }

    if (program == FC_PROGRAM_BAMSTATS && args->threads > 1) {
        pool.pool = hts_tpool_init(args->threads);
        if (pool.pool) {
            FC_LOG_INFO("Using %d threads", args->threads);
        } else {
            FC_LOG_WARN("Failed to create thread pool; continuing single-threaded.");
        }
    }

    // setup file/read iterator. This handles details of iterating over multiple files, directories, and/or file-of-filenames (fofn),
    bool region_mode = (program == FC_PROGRAM_BAMSTATS) && (args->region != NULL || args->bed != NULL);
    fc_iter_opts iter_opts = {
        .min_length = args->min_length,
        .max_length = args->max_length,
        .min_qscore = args->min_qscore,
        .dust = args->dust,
        .dust_t = (int)args->dust_t,
        .dust_w = (int)args->dust_w,
        .max_dust = (float)args->max_dust,
        .read_group = args->read_group,
        .tag_name = {args->tag_name[0], args->tag_name[1]},
        .tag_value = args->tag_value,
        .force_recalc_qual = args->force_recalc_qual,
        .region = region_mode ? args->region : NULL,
        .bed = region_mode ? args->bed : NULL,
        .thread_pool = (program == FC_PROGRAM_BAMSTATS && pool.pool) ? &pool : NULL
    };
    const char *bam_first_path = NULL;
    if (program == FC_PROGRAM_FASTCAT) {
        if (fc_iter_init(
                &file_iter,
                &args->input_paths,
                args->recurse,
                args->verbose,
                args->input_kind,
                args->fofn) != 0) {
            FC_LOG_ERROR("Failed to initialise file iterator.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
    } else {
        // BAMSTATS only takes a single BAM input, so we create a synthetic path list with just the one BAM
        // then bootstrap the iterator on that. This allows us to reuse the same iteration logic for both programs.
        char *bam_path_storage[] = {(char *)args->bam};
        fc_path_list input_paths = {
            .paths = bam_path_storage,
            .n_paths = 1
        };
        if (fc_iter_init(&file_iter, &input_paths, 0, args->verbose, FC_FILE_BAM, false) != 0
            || file_iter.n_files == 0) {
            FC_LOG_ERROR("No BAM inputs found.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (!fc_iter_next_file(&file_iter, &bam_first_path)) {
            FC_LOG_ERROR("No BAM inputs found.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (fc_iter_open_record_cursor(&file_iter, &iter_opts) != 0) {
            FC_LOG_ERROR("Failed to open BAM '%s'.", bam_first_path);
            status = EXIT_FAILURE;
            goto cleanup;
        }
    }

    // the write object controls dispatching reads to the various output sinks,
    // it tracks some knowledge of stats for the various sinks to use rather than
    // each sink having to calculate everything independently
    fc_stats_engine_config stats_cfg = {
        .track_alignment = (program == FC_PROGRAM_BAMSTATS),
        .track_polya = (program == FC_PROGRAM_BAMSTATS) ? args->poly_a : false,
        .polya_cover = (program == FC_PROGRAM_BAMSTATS) ? args->poly_a_cover : 0.0f,
        .polya_qual = (program == FC_PROGRAM_BAMSTATS) ? args->poly_a_qual : 0.0f,
        .polya_rev = (program == FC_PROGRAM_BAMSTATS) ? args->poly_a_rev : false
    };
    writer = fc_writer_create(args->demultiplex, &stats_cfg);
    if (!writer) {
        status = EXIT_FAILURE;
        goto cleanup;
    }

    // per-input-file summaries (e.g. FASTQ header stats, BAM flagstats)
    fc_writer_add_sink(writer, fc_perfile_sink_create(args->output_dir, args->sample));

    // histograms (read lengths, quality, etc.)
    fc_hist_sink_opts hist_opts = {
        .output_root = args->output_dir,
        .use_barcodes = args->demultiplex,
        .unmapped_override = (program == FC_PROGRAM_FASTCAT),
        .write_alignment_metrics = (program == FC_PROGRAM_BAMSTATS),
        .write_polya = (program == FC_PROGRAM_BAMSTATS) ? args->poly_a : false,
        .engines = writer->stats_engines,
        .engine_count = writer->stats_engine_count
    };
    fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));

    // counters (overall and per-barcode read counts, etc.)
    fc_counter_sink_opts counter_opts = {
        .output_root = args->output_dir,
        .sample = args->sample,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count
    };
    fc_writer_add_sink(writer, fc_counter_sink_create(&counter_opts));

    // alignment specific flags/statistics (BAMSTATS only)
    if (program == FC_PROGRAM_BAMSTATS) {
        fc_flagstats_sink_opts flag_opts = {
            .region_mode = region_mode,
            .output_root = args->output_dir,
            .sample = args->sample,
            .hdr = file_iter.bmi.hdr,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_flagstats_sink_create(&flag_opts));
    }

    // alignment coverage (BAMSTATS only, and only if --coverage specified)
    if (program == FC_PROGRAM_BAMSTATS && args->coverage) {        
        fc_coverage_sink_opts cov_opts = {
            .output_root = args->output_dir,
            .per_base = true,
            .use_cigar = false,
            .exclude_flags = -1,
            .include_flags = -1,
            .by_strand = true,
            .hdr = file_iter.bmi.hdr,
            .pool = pool.pool ? &pool : NULL,
            .bed_files = args->coverage_beds,
            .bed_names = args->coverage_names,
            .n_beds = args->n_coverage_beds,
            .thresholds = args->coverage_thresholds,
            .n_thresholds = args->n_coverage_thresholds,
            .segments = args->segments,
            .n_segments = args->n_segments,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_coverage_sink_create(&cov_opts));
    }

    // final output (FASTQ files for FASTCAT, and optionally BAM files if --write-bam specified)
    if (program == FC_PROGRAM_FASTCAT) {
        fc_output_sink_opts output_opts = {
            .output_root = args->demultiplex ? args->output_dir : NULL,
            .reads_per_file = args->reads_per_file,
            .write_bam = args->write_bam,
            .threads = args->threads,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_output_sink_create(&output_opts));
    }

    // line-by-line (read-by-read) detailed output.
    // This is optional since it can be very verbose.
    if (args->per_read) {
        fc_perread_sink_opts perread_opts = {
            .output_root = args->output_dir,
            .sample = args->sample,
            .include_alignment = (program == FC_PROGRAM_BAMSTATS),
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    }

    // global stats tracking across all files, for final summary output at the end
    fc_file_stats global_stats = {0};
    fc_file_stats_init(&global_stats);
    size_t nfile = file_iter.n_files;

    // go go gadget read iterator
    if (program == FC_PROGRAM_BAMSTATS) {
        status = process_current_file(
            &file_iter, writer, bam_first_path, file_iter.bmi.hdr, &iter_opts, true, false, &hard_error,
            &global_stats);
    } else {
        const char *path = NULL;
        while (fc_iter_next_file(&file_iter, &path)) {
            int file_status = process_current_file(
                &file_iter, writer, path, NULL, &iter_opts, false,
                (args->input_kind == FC_FILE_BAM), &hard_error, &global_stats);
            if (file_status != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
            }
            if (hard_error) {
                break;
            }
        }
    }

    fc_writer_finish(writer);

    uint64_t total_records =
        global_stats.failures[R_RECORD_OK]
        + global_stats.failures[R_TOO_LONG]
        + global_stats.failures[R_TOO_SHORT]
        + global_stats.failures[R_LOW_QUALITY]
        + global_stats.failures[R_DUST_MASKED];
    FC_LOG_INFO("Processed %" PRIu64 " records in %zu files.", total_records, nfile);
    if (status != EXIT_SUCCESS) {
        FC_LOG_WARN("Error processing files.");
    }
    if (!args->force_error && !hard_error) {
        status = EXIT_SUCCESS;
    }

    FC_LOG_INFO("Parsing/filtering summary:");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        FC_LOG_INFO("%s\t%" PRIu64, failure_type[i], global_stats.failures[i]);
    }
    clock_t end = clock();
    FC_LOG_INFO("Total CPU time: %fs", (double)(end - begin) / CLOCKS_PER_SEC);

cleanup:
    fc_writer_destroy(writer);
    fc_iter_close(&file_iter);
    if (pool.pool) {
        hts_tpool_destroy(pool.pool);
    }
    fc_log_destroy();
    return status;
}
