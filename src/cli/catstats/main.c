#include "cli/catstats/main.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "htslib/thread_pool.h"

#include "core/fc_bam_utils.h"
#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_iter.h"
#include "core/fc_log.h"
#include "core/fc_report.h"
#include "core/fc_report_arg_table.h"
#include "core/fc_sinks.h"
#include "core/fc_writer.h"

static bool filters_requested(const fc_catstats_args *args)
{
    return args->min_length != 0 || args->max_length != (size_t)-1 || args->min_qscore != 0.0f
        || args->dust || args->read_group != NULL || args->tag_name[0] != '\0';
}

static void log_bam_input_error(fc_bam_input_status status, const char *path)
{
    switch (status) {
    case FC_BAM_INPUT_OPEN_FAILED:
        FC_LOG_ERROR("Unable to open BAM '%s' during header preflight.", path);
        break;
    case FC_BAM_INPUT_HEADER_FAILED:
        FC_LOG_ERROR("Unable to read BAM header from '%s'.", path);
        break;
    case FC_BAM_INPUT_REFERENCE_MISMATCH:
        FC_LOG_ERROR("Reference dictionary in '%s' does not match earlier BAM input.", path);
        break;
    case FC_BAM_INPUT_HEADER_COPY_FAILED:
        FC_LOG_ERROR("Unable to retain reference header from '%s'.", path);
        break;
    case FC_BAM_INPUT_OK:
        break;
    }
}

static fc_input_resolve_status resolve_program_inputs(
    fc_program_kind program,
    const fc_catstats_args *args,
    fc_resolved_inputs *inputs)
{
    bool filter_kind = true;
    fc_file_kind requested_kind = FC_FILE_FASTQ;
    if (program == FC_PROGRAM_BAMSTATS) {
        requested_kind = FC_FILE_BAM;
    } else if (program == FC_PROGRAM_QC) {
        filter_kind = args->input_kind_explicit;
        requested_kind = args->input_kind;
    }
    fc_input_resolve_status status = fc_resolve_input_paths(
        &args->input_paths,
        args->recurse,
        args->fofn,
        filter_kind,
        requested_kind,
        inputs);
    if (status != FC_INPUT_RESOLVE_OK) {
        return status;
    }
    if (program == FC_PROGRAM_BAMSTATS && inputs->paths.n_paths != 1) {
        return FC_INPUT_RESOLVE_INVALID_ARGUMENT;
    }
    if (filter_kind) {
        return FC_INPUT_RESOLVE_OK;
    }
    if (inputs->found_kinds == FC_FILE_KINDS_FASTQ) {
        inputs->kind = FC_FILE_FASTQ;
        return FC_INPUT_RESOLVE_OK;
    }
    if (inputs->found_kinds == FC_FILE_KINDS_BAM) {
        inputs->kind = FC_FILE_BAM;
        return FC_INPUT_RESOLVE_OK;
    }
    return FC_INPUT_RESOLVE_MIXED_FORMAT;
}

static int process_current_file(
    fc_iter *file_iter,
    fc_writer *writer,
    const char *path,
    const sam_hdr_t *hdr_for_records,
    const fc_iter_opts *iter_opts,
    bool enforce_coordinate_order,
    bool *hard_error,
    fc_file_stats *global_stats)
{
    int status = EXIT_SUCCESS;

    fc_writer_begin_file(writer, path);

    if (!file_iter->cursor_open && fc_iter_open_record_cursor(file_iter, iter_opts) != 0) {
        status = EXIT_FAILURE;
        fc_iter_merge_file_stats(file_iter, global_stats);
        fc_writer_end_file(writer, path, fc_iter_stats(file_iter));
        return status;
    }

    const sam_hdr_t *record_header = hdr_for_records;
    if (record_header == NULL && file_iter->expected_kind == FC_FILE_BAM) {
        record_header = file_iter->bmi.hdr;
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

    fc_bam_coordinate_order coordinate_order = { 0 };
    fc_iter_next_status next_rc = FC_ITER_RECORD;
    while ((next_rc = fc_iter_next_record(file_iter, fr)) == FC_ITER_RECORD) {
        if (enforce_coordinate_order
            && !fc_bam_coordinate_order_check(&coordinate_order, fr->bam)) {
            FC_LOG_ERROR("Aligned BAM stream is not coordinate sorted.");
            if (hard_error) {
                *hard_error = true;
            }
            status = EXIT_FAILURE;
            break;
        }
        fc_writer_add_read(writer, fr, path, record_header);
    }
    if (next_rc == FC_ITER_ERROR) {
        status = EXIT_FAILURE;
        if (fc_iter_has_hard_error(file_iter) && hard_error != NULL) {
            *hard_error = true;
        }
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

    const char *log_module = "FASTCAT";
    if (program == FC_PROGRAM_BAMSTATS) {
        log_module = "BAMSTATS";
    } else if (program == FC_PROGRAM_QC) {
        log_module = "FASTCAT-QC";
    }
    fc_log_config log_cfg = {
        .module = log_module,
        .level = args->log_level,
        .target = args->log_file,
        .color = args->log_mono ? FC_LOG_COLOR_OFF : FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);
    if (program == FC_PROGRAM_FASTCAT || program == FC_PROGRAM_BAMSTATS) {
        FC_LOG_WARN("'fastcat %s' is retained for compatibility; prefer 'fastcat qc' for new QC reports.",
            program == FC_PROGRAM_BAMSTATS ? "bam" : "fastq");
    }

    // initialised below if we successfully create outputs
    htsThreadPool pool = { NULL, 0 };
    fc_writer *writer = NULL;
    fc_alignment_summary *alignment_summary = NULL;
    fc_coverage_summaries *coverage_summaries = NULL;
    fc_iter file_iter = { 0 };
    fc_resolved_inputs resolved_inputs = { 0 };
    fc_bam_input_info bam_input_info = { 0 };
    const sam_hdr_t *reference_header = NULL;
    bool hard_error = false;
    bool bam_mode = false;
    bool has_reference_bam = false;
    bool calculate_coverage = false;
    bool direct_reads_output = program == FC_PROGRAM_QC && args->reads_output != NULL
        && !args->demultiplex;
    bool output_created = false;
    bool output_published = false;
    char *final_output_dir = NULL;
    char *temporary_output_dir = NULL;
    char *temporary_reads_output = NULL;
    const char *output_dir = NULL;
    fc_input_resolve_status resolve_status = resolve_program_inputs(
        program, args, &resolved_inputs);
    if (resolve_status != FC_INPUT_RESOLVE_OK) {
        if (program == FC_PROGRAM_BAMSTATS
            && resolve_status == FC_INPUT_RESOLVE_INVALID_ARGUMENT
            && resolved_inputs.paths.n_paths > 1) {
            FC_LOG_ERROR("fastcat bam accepts exactly one BAM input.");
        } else if (program == FC_PROGRAM_FASTCAT
            && resolve_status == FC_INPUT_RESOLVE_EMPTY
            && resolved_inputs.found_kinds == FC_FILE_KINDS_BAM) {
            FC_LOG_ERROR("No FASTQ input files found; use 'fastcat qc' for BAM input.");
        } else {
            FC_LOG_ERROR("Failed to resolve input paths: %s.",
                fc_input_resolve_status_string(resolve_status));
        }
        status = EXIT_FAILURE;
        goto cleanup;
    }
    fc_file_kind input_kind = resolved_inputs.kind;
    bam_mode = input_kind == FC_FILE_BAM;
    if (program == FC_PROGRAM_QC && !bam_mode
        && (args->region || args->bed || args->coverage || args->coverage_configured || args->poly_a)) {
        FC_LOG_ERROR("--region, --bed, coverage, and poly-A options require BAM input.");
        status = EXIT_FAILURE;
        goto cleanup;
    }
    size_t final_output_len = strlen(args->output_dir);
    while (final_output_len > 1 && args->output_dir[final_output_len - 1] == '/') {
        --final_output_len;
    }
    final_output_dir = fc_strndup(args->output_dir, final_output_len, "output directory");

#ifdef NOTHREADS
    if (args->threads != 1) {
        fprintf(stderr, "--threads set to %d, but threading not supported by this build.\n", args->threads);
    }
#endif

    if (final_output_dir[0] == '\0') {
        FC_LOG_ERROR("Output directory must not be empty.");
        status = EXIT_FAILURE;
        goto cleanup;
    }
    temporary_output_dir = fc_asprintf(
        "temporary output directory", "%s.tmp", final_output_dir);
    output_dir = temporary_output_dir;
    if (file_exists(final_output_dir)) {
        FC_LOG_ERROR("Output directory '%s' already exists.", final_output_dir);
        status = EXIT_FAILURE;
        goto cleanup;
    }
    if (file_exists(temporary_output_dir)) {
        FC_LOG_ERROR("Temporary output directory '%s' already exists.", temporary_output_dir);
        status = EXIT_FAILURE;
        goto cleanup;
    }
    if (mkdir_hier(temporary_output_dir, true) == -1) {
        FC_LOG_ERROR("Unable to create temporary output directory '%s'.", temporary_output_dir);
        status = EXIT_FAILURE;
        goto cleanup;
    }
    output_created = true;
    if (direct_reads_output) {
        const char *reads_extension = args->write_bam ? ".bam" : ".fastq.gz";
        char *reads_filename = fc_asprintf(
            "read output filename", "%s%s", args->reads_output, reads_extension);
        temporary_reads_output = fc_join_path(output_dir, reads_filename);
        free(reads_filename);
        if (ensure_parent_dir_exists(temporary_reads_output) != 0) {
            FC_LOG_ERROR("Unable to create parent directory for read output prefix '%s'.",
                args->reads_output);
            status = EXIT_FAILURE;
            goto cleanup;
        }
    }

    if (bam_mode && args->threads > 1) {
        pool.pool = hts_tpool_init(args->threads);
        if (pool.pool) {
            FC_LOG_INFO("Using %d threads", args->threads);
        } else {
            FC_LOG_WARN("Failed to create thread pool; continuing single-threaded.");
        }
    }

    // setup file/read iterator. This handles details of iterating over multiple files, directories, and/or file-of-filenames (fofn),
    bool region_mode = bam_mode && (args->region != NULL || args->bed != NULL);
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
        .thread_pool = (bam_mode && pool.pool) ? &pool : NULL,
        .use_barcodes = args->demultiplex
    };
    if (fc_iter_init(&file_iter, &resolved_inputs, args->verbose) != 0) {
        FC_LOG_ERROR("Failed to initialise file iterator.");
        status = EXIT_FAILURE;
        goto cleanup;
    }
    if (bam_mode) {
        const char *failed_path = NULL;
        fc_bam_input_status preflight = fc_bam_input_preflight(
            (const char *const *)file_iter.files, file_iter.n_files, &bam_input_info, &failed_path);
        if (preflight != FC_BAM_INPUT_OK) {
            log_bam_input_error(preflight, failed_path);
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (bam_input_info.has_stdin && file_iter.n_files != 1) {
            FC_LOG_ERROR("BAM stdin cannot be combined with file inputs.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (program == FC_PROGRAM_QC && (args->coverage || args->coverage_configured)
            && file_iter.n_files != 1) {
            FC_LOG_ERROR("Coverage is supported only with one BAM input.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (bam_input_info.has_stdin) {
            if (!fc_iter_next_file(&file_iter, NULL)
                || fc_iter_open_record_cursor(&file_iter, &iter_opts) != 0) {
                FC_LOG_ERROR("Failed to open BAM stream.");
                status = EXIT_FAILURE;
                goto cleanup;
            }
            fc_bam_input_status header_status = fc_bam_input_info_add_header(
                &bam_input_info, file_iter.bmi.hdr);
            if (header_status != FC_BAM_INPUT_OK) {
                log_bam_input_error(header_status, file_iter.current);
                status = EXIT_FAILURE;
                goto cleanup;
            }
        }
        reference_header = bam_input_info.reference_header;
        has_reference_bam = bam_input_info.has_references;
        if (program == FC_PROGRAM_QC && has_reference_bam && direct_reads_output) {
            FC_LOG_ERROR("--reads-output is not supported for reference-bearing BAM QC input.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (program == FC_PROGRAM_QC && has_reference_bam && filters_requested(args)) {
            FC_LOG_WARN("Read filtering is disabled for reference-bearing BAM QC input.");
            fc_iter_opts_disable_filters(&iter_opts);
            if (file_iter.cursor_open) {
                fc_iter_opts_disable_filters(&file_iter.bmi.opts);
            }
        }
        if (!has_reference_bam && (args->coverage || args->coverage_configured)) {
            FC_LOG_ERROR("Coverage requires a reference-bearing BAM input.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (region_mode && !has_reference_bam) {
            FC_LOG_ERROR("--region and --bed require a reference-bearing BAM input.");
            status = EXIT_FAILURE;
            goto cleanup;
        }
        if (has_reference_bam && file_iter.n_files == 1
            && fc_bam_header_is_coordinate_sorted(reference_header)) {
            calculate_coverage = true;
        } else if (args->coverage || args->coverage_configured) {
            FC_LOG_ERROR("Coverage requires a coordinate-sorted BAM input.");
            status = EXIT_FAILURE;
            goto cleanup;
        } else if (has_reference_bam) {
            FC_LOG_WARN("Coverage summary omitted because BAM input is not coordinate sorted.");
        }
        if (program == FC_PROGRAM_QC && has_reference_bam) {
            FC_LOG_WARN("Read output is disabled for reference-bearing BAM QC input.");
        }
    }

    // the write object controls dispatching reads to the various output sinks,
    // it tracks some knowledge of stats for the various sinks to use rather than
    // each sink having to calculate everything independently
    fc_stats_engine_config stats_cfg = {
        .track_alignment = bam_mode,
        .track_polya = bam_mode ? args->poly_a : false,
        .polya_cover = bam_mode ? args->poly_a_cover : 0.0f,
        .polya_qual = bam_mode ? args->poly_a_qual : 0.0f,
        .polya_rev = bam_mode ? args->poly_a_rev : false
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
    fc_accuracy_qscore_2d_config accuracy_qscore_cfg = {
        .kll_k = args->length_qscore_k,
        .qgrid_points = args->length_qscore_qgrid,
        .smoothing_points = args->length_qscore_qsmooth
    };
    fc_time_bin_config time_bin_cfg = {
        .bin_size_s = args->time_bin_size_s,
        .kll_k = args->time_bin_k
    };
    bool bam_input = bam_mode;
    if (bam_input) {
        size_t n_refs = reference_header != NULL
            ? (size_t)reference_header->n_targets
            : 0;
        size_t barcode_count = args->demultiplex ? FC_MAX_BARCODES : 1;
        alignment_summary = fc_alignment_summary_create(
            n_refs,
            region_mode,
            args->demultiplex,
            barcode_count);
    }
    fc_writer_config writer_cfg = {
        .use_barcodes = args->demultiplex,
        .alignment_summary = alignment_summary,
        .stats = &stats_cfg,
        .run_summary = true,
        .channel_summary = true,
        .length_qscore = &length_qscore_cfg,
        .speed_qscore = &speed_qscore_cfg,
        .accuracy_qscore = bam_input ? &accuracy_qscore_cfg : NULL,
        .time_bin = &time_bin_cfg
    };
    writer = fc_writer_create(&writer_cfg);
    if (!writer) {
        status = EXIT_FAILURE;
        goto cleanup;
    }
    iter_opts.run_summary = writer->run_summary;
    if (file_iter.cursor_open) {
        file_iter.bmi.opts.run_summary = writer->run_summary;
    }

    bool write_legacy_outputs = args->legacy_outputs;
    if (write_legacy_outputs) {
        // Per-input-file summaries (e.g. FASTQ header stats, BAM flagstats).
        fc_writer_add_sink(writer, fc_perfile_sink_create(output_dir, args->sample));

        // histograms (read lengths, quality, etc.)
        fc_hist_sink_opts hist_opts = {
            .output_root = output_dir,
            .use_barcodes = args->demultiplex,
            .write_alignment_metrics = bam_mode,
            .write_polya = bam_mode ? args->poly_a : false,
            .engines = writer->stats_engines,
            .engine_count = writer->stats_engine_count
        };
        fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));

        // counters (overall and per-barcode read counts, etc.)
        fc_counter_sink_opts counter_opts = {
            .output_root = output_dir,
            .sample = args->sample,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_counter_sink_create(&counter_opts));
    }

    // Coverage snapshots are consumed by later sink finalisers, so register their producer first.
    if (calculate_coverage) {
        coverage_summaries = fc_coverage_summaries_create(writer->barcode_count);
        fc_coverage_sink_opts cov_opts = {
            .output_root = output_dir,
            .write_traces = args->coverage,
            .write_tables = program == FC_PROGRAM_BAMSTATS
                ? args->coverage && write_legacy_outputs
                : write_legacy_outputs,
            .use_cigar = false,
            .exclude_flags = -1,
            .include_flags = -1,
            .by_strand = true,
            .hdr = reference_header,
            .pool = pool.pool ? &pool : NULL,
            .bed_files = args->coverage_beds,
            .bed_names = args->coverage_names,
            .n_beds = args->n_coverage_beds,
            .thresholds = args->coverage_thresholds,
            .n_thresholds = args->n_coverage_thresholds,
            .segments = args->segments,
            .n_segments = args->n_segments,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count,
            .summaries = coverage_summaries
        };
        fc_writer_add_sink(writer, fc_coverage_sink_create(&cov_opts));
    }

    fc_summary_json_sink_opts summary_json_opts = {
        .output_root = output_dir,
        .sample = args->sample,
        .use_barcodes = args->demultiplex,
        .barcode_count = writer->barcode_count,
        .bam_input = bam_input,
        .region_mode = region_mode,
        .region_label = region_mode ? (args->region ? args->region : "bed_regions") : NULL,
        .hdr = reference_header,
        .engines = writer->stats_engines,
        .engine_count = writer->stats_engine_count,
        .run_summary = writer->run_summary,
        .channel_summary = writer->channel_summary,
        .length_qscore = writer->length_qscore,
        .speed_qscore = writer->speed_qscore,
        .accuracy_qscore = writer->accuracy_qscore,
        .time_bin = writer->time_bin,
        .alignment_summary = alignment_summary,
        .coverage_summaries = coverage_summaries
    };
    fc_writer_add_sink(writer, fc_summary_json_sink_create(&summary_json_opts));

    // Alignment-specific flags/statistics for BAM input.
    if (bam_mode && write_legacy_outputs) {
        const char *region_label = NULL;
        if (region_mode) {
            region_label = args->region ? args->region : "bed_regions";
        }
        fc_flagstats_sink_opts flag_opts = {
            .region_mode = region_mode,
            .region_label = region_label,
            .output_root = output_dir,
            .sample = args->sample,
            .hdr = reference_header,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count,
            .alignment_summary = alignment_summary
        };
        fc_writer_add_sink(writer, fc_flagstats_sink_create(&flag_opts));
    }

    // Reference-bearing BAM is analysis-only: concatenating selected coordinate
    // regions or multiple BAMs would not produce a valid alignment stream.
    bool emit_reads = program != FC_PROGRAM_BAMSTATS && !has_reference_bam;
    if (emit_reads) {
        fc_output_sink_opts output_opts = {
            .output_root = args->demultiplex ? output_dir : NULL,
            .stream_path = direct_reads_output ? temporary_reads_output : NULL,
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
            .output_root = output_dir,
            .sample = args->sample,
            .include_alignment = bam_mode,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    }

    // global stats tracking across all files, for final summary output at the end
    fc_file_stats global_stats = { 0 };
    fc_file_stats_init(&global_stats);
    size_t nfile = file_iter.n_files;

    const char *path = NULL;
    if (file_iter.cursor_open) {
        status = process_current_file(
            &file_iter, writer, file_iter.current, reference_header, &iter_opts,
            bam_mode && has_reference_bam, &hard_error, &global_stats);
    }
    while (!hard_error && fc_iter_next_file(&file_iter, &path)) {
        int file_status = process_current_file(
            &file_iter, writer, path, NULL, &iter_opts, false, &hard_error, &global_stats);
        if (file_status != EXIT_SUCCESS) {
            status = EXIT_FAILURE;
        }
    }

    // Match fastq-mode semantics: seeing records without quality marks the run
    // as failed unless caller opts out via the final non-force-error override.
    if (global_stats.failures[F_QUAL_MISSING] > 0) {
        status = EXIT_FAILURE;
    }

    fc_writer_finish(writer);

    if (args->html_report) {
        fc_report_arg_table report_arg_table = fc_report_arg_table_from_catstats(
            program, args, input_kind, nfile);
        fc_report_opts report_opts = {
            .output_root = output_dir,
            .sample = args->sample,
            .run_config = report_arg_table.entries,
            .run_config_count = report_arg_table.n_entries,
            .bam_input = bam_input,
            .use_barcodes = args->demultiplex,
            .barcode_count = writer->barcode_count,
            .engines = writer->stats_engines,
            .engine_count = writer->stats_engine_count,
            .run_summary = writer->run_summary,
            .channel_summary = writer->channel_summary,
            .length_qscore = writer->length_qscore,
            .speed_qscore = writer->speed_qscore,
            .accuracy_qscore = writer->accuracy_qscore,
            .time_bin = writer->time_bin,
            .alignment_summary = alignment_summary,
            .coverage_summaries = coverage_summaries
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
    fc_resolved_inputs_free(&resolved_inputs);
    fc_writer_destroy(writer);
    fc_alignment_summary_destroy(alignment_summary);
    fc_coverage_summaries_destroy(coverage_summaries);
    fc_iter_close(&file_iter);
    fc_bam_input_info_destroy(&bam_input_info);
    if (pool.pool) {
        hts_tpool_destroy(pool.pool);
    }
    if (output_created && status == EXIT_SUCCESS) {
        if (rename(temporary_output_dir, final_output_dir) != 0) {
            FC_LOG_ERROR("Unable to publish output directory '%s': %s.",
                final_output_dir, strerror(errno));
            status = EXIT_FAILURE;
        } else {
            output_published = true;
        }
    }
    if (output_created && !output_published) {
        if (args->keep_failed_output) {
            FC_LOG_WARN("Keeping partial output in '%s'.", temporary_output_dir);
        } else if (fc_remove_tree(temporary_output_dir) != 0) {
            FC_LOG_ERROR("Unable to remove partial output directory '%s': %s.",
                temporary_output_dir, strerror(errno));
        }
    }
    free(temporary_output_dir);
    free(temporary_reads_output);
    free(final_output_dir);
    fc_log_destroy();
    return status;
}
