#include "cli/catstats/main.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "htslib/thread_pool.h"

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_iter.h"
#include "core/fc_log.h"
#include "core/fc_report.h"
#include "core/fc_report_arg_table.h"
#include "core/fc_sinks.h"
#include "core/fc_writer.h"

static int process_current_file(
    fc_iter *file_iter,
    fc_writer *writer,
    const char *path,
    const sam_hdr_t *hdr_for_records,
    const fc_iter_opts *iter_opts,
    bool cursor_already_open,
    bool enforce_unmapped,
    bool *hard_error,
    fc_file_stats *global_stats)
{
    int status = EXIT_SUCCESS;

    fc_writer_begin_file(writer, path);

    if (!cursor_already_open && fc_iter_open_record_cursor(file_iter, iter_opts) != 0) {
        status = EXIT_FAILURE;
        fc_iter_merge_file_stats(file_iter, global_stats);
        fc_writer_end_file(writer, path, fc_iter_stats(file_iter));
        return status;
    }

    fc_read *fr = fc_read_init();
    if (!fr) {
        FC_LOG_ERROR("Failed to allocate read buffer.");
        if (hard_error) {
            *hard_error = true;
        }
        status = EXIT_FAILURE;
        goto finish_file;
    }

    fc_iter_next_status next_rc = FC_ITER_RECORD;
    while ((next_rc = fc_iter_next_record(file_iter, fr)) == FC_ITER_RECORD) {
        if (enforce_unmapped && !fc_is_unmapped(fr->bam)) {
            FC_LOG_ERROR(
                "Input BAM contains mapped reads (e.g. '%s'). "
                "fastcat fastq with --input-fmt bam only accepts unaligned BAM input.",
                bam_get_qname(fr->bam));
            if (hard_error) {
                *hard_error = true;
            }
            status = EXIT_FAILURE;
            break;
        }
        fc_writer_add_read(writer, fr, path, hdr_for_records);
    }
    if (next_rc == FC_ITER_ERROR) {
        status = EXIT_FAILURE;
    }

    fc_read_destroy(fr);

finish_file:
    if (fc_iter_finish_file(file_iter) != 0) {
        status = EXIT_FAILURE;
    }
    fc_iter_merge_file_stats(file_iter, global_stats);
    fc_writer_end_file(writer, path, fc_iter_stats(file_iter));
    return status;
}

int fc_run_program(fc_program_kind program, const fc_catstats_args *args)
{
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
    htsThreadPool pool = { NULL, 0 };
    fc_writer *writer = NULL;
    fc_iter file_iter = { 0 };
    bool hard_error = false;

#ifdef NOTHREADS
    if (args->threads != 1) {
        fprintf(stderr, "--threads set to %d, but threading not supported by this build.\n", args->threads);
    }
#endif

    if (mkdir_hier(args->output_dir, true) == -1) {
        FC_LOG_ERROR("Output directory '%s' already exists or cannot be created.", args->output_dir);
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
        .tag_name = { args->tag_name[0], args->tag_name[1] },
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
                args->fofn)
            != 0) {
            FC_LOG_ERROR("Failed to initialise file iterator.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (file_iter.n_files == 0) {
            FC_LOG_ERROR("No usable input files found.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
    } else {
        // BAMSTATS only takes a single BAM input, so we create a synthetic path list with just the one BAM
        // then bootstrap the iterator on that. This allows us to reuse the same iteration logic for both programs.
        char *bam_path_storage[] = { (char *)args->bam };
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
    fc_length_qscore_2d_config length_qscore_cfg = {
        .kll_k = args->length_qscore_k,
        .qgrid_points = args->length_qscore_qgrid,
        .smoothing_points = args->length_qscore_qsmooth
    };
    fc_speed_qscore_2d_config speed_qscore_cfg = {
        .kll_k = args->length_qscore_k,
        .qgrid_points = args->length_qscore_qgrid,
        .smoothing_points = args->length_qscore_qsmooth
    };
    fc_time_bin_config time_bin_cfg = {
        .bin_size_s = args->time_bin_size_s,
        .kll_k = args->time_bin_k
    };
    fc_writer_config writer_cfg = {
        .use_barcodes = args->demultiplex,
        .stats = &stats_cfg,
        .run_summary = true,
        .channel_summary = true,
        .length_qscore = &length_qscore_cfg,
        .speed_qscore = &speed_qscore_cfg,
        .time_bin = &time_bin_cfg
    };
    writer = fc_writer_create(&writer_cfg);
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

    fc_run_summary_sink_opts run_summary_opts = {
        .output_root = args->output_dir,
        .sample = args->sample,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .bam_input = (program == FC_PROGRAM_BAMSTATS) || (args->input_kind == FC_FILE_BAM),
        .unmapped_override = (program == FC_PROGRAM_FASTCAT),
        .engines = writer->stats_engines,
        .engine_count = writer->stats_engine_count,
        .model = writer->run_summary
    };
    fc_writer_add_sink(writer, fc_run_summary_sink_create(&run_summary_opts));

    fc_channel_summary_sink_opts channel_opts = {
        .output_root = args->output_dir,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .model = writer->channel_summary
    };
    fc_writer_add_sink(writer, fc_channel_summary_sink_create(&channel_opts));

    // read-length vs qscore 2D density table derived from writer-owned KLL sketches
    fc_length_qscore_2d_sink_opts lq2d_opts = {
        .output_root = args->output_dir,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .model = writer->length_qscore
    };
    fc_writer_add_sink(writer, fc_length_qscore_2d_sink_create(&lq2d_opts));

    fc_speed_qscore_2d_sink_opts qts2d_opts = {
        .output_root = args->output_dir,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .model = writer->speed_qscore
    };
    fc_writer_add_sink(writer, fc_speed_qscore_2d_sink_create(&qts2d_opts));

    // time-binned summaries for yield and boxplot-style temporal QC
    fc_time_bin_sink_opts time_bin_opts = {
        .output_root = args->output_dir,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .model = writer->time_bin
    };
    fc_writer_add_sink(writer, fc_time_bin_sink_create(&time_bin_opts));

    // alignment specific flags/statistics (BAMSTATS only)
    if (program == FC_PROGRAM_BAMSTATS) {
        const char *region_label = NULL;
        if (region_mode) {
            region_label = args->region ? args->region : "bed_regions";
        }
        fc_flagstats_sink_opts flag_opts = {
            .region_mode = region_mode,
            .region_label = region_label,
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
    fc_file_stats global_stats = { 0 };
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
                (args->input_kind != FC_FILE_FASTQ), &hard_error, &global_stats);
            if (file_status != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
            }
            if (hard_error) {
                break;
            }
        }
    }

    // Match fastq-mode semantics: seeing records without quality marks the run
    // as failed unless caller opts out via the final non-force-error override.
    if (global_stats.failures[F_QUAL_MISSING] > 0) {
        status = EXIT_FAILURE;
    }

    fc_writer_finish(writer);

    if (args->html_report) {
        fc_report_arg_table report_arg_table = fc_report_arg_table_from_catstats(program, args, nfile);
        fc_report_opts report_opts = {
            .output_root = args->output_dir,
            .sample = args->sample,
            .run_config = report_arg_table.entries,
            .run_config_count = report_arg_table.n_entries,
            .bam_input = (program == FC_PROGRAM_BAMSTATS) || (args->input_kind == FC_FILE_BAM),
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count,
            .unmapped_override = (program == FC_PROGRAM_FASTCAT),
            .engines = writer->stats_engines,
            .engine_count = writer->stats_engine_count,
            .run_summary = writer->run_summary,
            .channel_summary = writer->channel_summary,
            .length_qscore = writer->length_qscore,
            .speed_qscore = writer->speed_qscore,
            .time_bin = writer->time_bin
        };
        if (fc_write_html_reports(&report_opts) != 0) {
            hard_error = true;
            status = EXIT_FAILURE;
        }
        fc_report_arg_table_destroy(&report_arg_table);
    }

    uint64_t seen_records = global_stats.failures[R_RECORD_SEEN];
    uint64_t used_records = global_stats.failures[R_RECORD_OK];
    FC_LOG_INFO(
        "Processed %" PRIu64 " records (%" PRIu64 " passed filters) in %zu files.",
        seen_records, used_records, nfile);
    if (status != EXIT_SUCCESS) {
        FC_LOG_WARN("Error processing files.");
    }
    // for fastq-mode parsing/stream errors we allow isolated broken files to
    // not bring everything to a halt, unless --force-error is given
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
