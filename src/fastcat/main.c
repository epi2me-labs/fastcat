#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdlib.h>

#include "../common.h"
#include "../fc_log.h"
#include "args.h"
#include "parsing.h"
#include "../fc_iter.h"
#include "../fc_sinks.h"
#include "../fc_writer.h"


const char filetypes[4][10] = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};
size_t nfiletypes = 4;


int main(int argc, char **argv) {
    arguments_t args = parse_arguments(argc, argv);
    fc_log_config log_cfg = {
        .module = "FASTCAT",
        .level = args.log_level,
        .target = args.log_file,
        .color = args.log_mono ? FC_LOG_COLOR_OFF : FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    int can_output = can_make_dir(args.output_dir);
    if (can_output != 1) {
        FC_LOG_ERROR("Output directory '%s' already exists or cannot be created.", args.output_dir);
        exit(EXIT_FAILURE);
    }
    if (mkdir_hier(args.output_dir) == -1) {
        FC_LOG_ERROR("Failed to create output directory '%s'.", args.output_dir);
        exit(EXIT_FAILURE);
    }

    fc_stats_engine_config stats_cfg = {
        .track_alignment = false,
        .track_polya = false,
        .polya_cover = 0.0f,
        .polya_qual = 0.0f,
        .polya_rev = false
    };

    // derive default paths under the output root
    fc_writer *writer = fc_writer_create(args.demultiplex, &stats_cfg);
    if (writer == NULL) exit(1);

    fc_output_sink_opts output_opts = {
        // when demultiplexing, output should go to output_dir; else to stdout
        .output_root = args.demultiplex ? args.output_dir : NULL,
        .reads_per_file = args.reads_per_file,
        .reheader = args.reheader,
        .write_bam = args.write_bam,
        .threads = args.threads,
        .barcode_count = writer->barcode_count
    };
    fc_writer_add_sink(writer, fc_output_sink_create(&output_opts));

    fc_hist_sink_opts hist_opts = {
        .output_root = args.output_dir,
        .use_barcodes = args.demultiplex,
        .unmapped_override = true,
        .write_alignment_metrics = false,
        .write_polya = false,
        .engines = writer->stats_engines,
        .engine_count = writer->stats_engine_count
    };
    fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));
    if (args.per_read) {
        fc_perread_sink_opts perread_opts = {
            .output_root = args.output_dir,
            .sample = args.sample,
            .include_alignment = false,
            .use_barcodes = args.demultiplex,
            .barcode_count = writer->barcode_count
        };
        fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    }
    fc_writer_add_sink(writer, fc_perfile_sink_create(args.output_dir, args.sample));
    fc_counter_sink_opts counter_opts = {
        .output_root = args.output_dir,
        .sample = args.sample,
        .use_barcodes = args.demultiplex,
        .barcode_count = writer->barcode_count
    };
    fc_writer_add_sink(writer, fc_counter_sink_create(&counter_opts));

    fc_iter_opts iter_opts = {
        .min_length = args.min_length,
        .max_length = args.max_length,
        .min_qscore = args.min_qscore,
        .dust = args.dust,
        .dust_t = args.dust_t,
        .dust_w = args.dust_w,
        .max_dust = args.max_dust,
        .read_group = NULL,
        .tag_name = { 0, 0 },
        .tag_value = 0,
        .force_recalc_qual = false,
        .region = NULL,
        .bed = NULL
    };
    fc_file_stats global_stats = {0};
    fc_file_stats_init(&global_stats);

    // inputs prepared by parse_arguments (stdin handled there when "-" is given)
    char **input_paths = args.files;
    size_t n_inputs = args.n_files;
    size_t nfile = 0;

    int status = 0;
    fc_iter file_iter = {0};
    if (fc_iter_init(&file_iter, input_paths, n_inputs, args.recurse, args.verbose) != 0) {
        FC_LOG_ERROR("Failed to initialise file iterator.");
        exit(EXIT_FAILURE);
    }
    nfile = file_iter.n_files;

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
            fc_writer_add_read(writer, fr, f->path, NULL);
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
    FC_LOG_INFO("Processed %" PRIu64 " records in %zu files.", total_records, nfile);
    if (status != EXIT_SUCCESS) {
        FC_LOG_WARN("Error processing files.");
    }
    if (!args.force_error) {
        status = EXIT_SUCCESS;
    }

    FC_LOG_INFO("Parsing/filtering summary:");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        FC_LOG_INFO("%s\t%" PRIu64, failure_type[i], global_stats.failures[i]);
    }
    fc_writer_destroy(writer);
    fc_log_destroy();
    destroy_args(&args);
    return status;
}
