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
#include "args.h"
#include "parsing.h"
#include "../fc_iter.h"
#include "../fc_sinks.h"
#include "../fc_writer.h"


const char filetypes[4][10] = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};
size_t nfiletypes = 4;

static char* strip_trailing_slash(const char* input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (len > 0 && input[len - 1] == '/') len--;
    char* output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
}

// defined below -- recursion
int process_file(char* fname, fc_writer* writer, const fc_fastq_iter_opts* iter_opts,
        uint64_t* failures, arguments_t *args, int recurse);

int process_dir(const char *name, fc_writer* writer, const fc_fastq_iter_opts* iter_opts,
        uint64_t* failures, arguments_t *args, int recurse) {
    int status = 0;
    DIR *dir;
    struct dirent *entry;
    char* search;

    // read all files in directory
    if (!(dir = opendir(name))) {
        fprintf(stderr, "ERROR  : could not process directory %s: %s\n", name, strerror(errno));
        return errno;
    }
    while ((entry = readdir(dir)) != NULL) {
        char *path = calloc(strlen(name) + strlen(entry->d_name) + 2, sizeof(char));
        sprintf(path, "%s/%s", name, entry->d_name);
        if ((entry->d_type == DT_DIR) && (recurse != 0)) {
            // skip
        } else {
            for (size_t i=0; i<nfiletypes; ++i) {
                search = strstr(entry->d_name, filetypes[i]);
                if (search != NULL) {
                    if (args->verbose) {
                        fprintf(stderr, "Processing %s\n", path);
                    }
                    int rtn = process_file(path, writer, iter_opts, failures, args, recurse - 1);
                    status = max(status, rtn);
                    break;
                }
            }
        }
        free(path);
    }
    closedir(dir);

    // start again and look at child directories
    if (!(dir = opendir(name))) {
        fprintf(stderr, "ERROR  : could not process directory %s: %s\n", name, strerror(errno));
        return errno;
    }
    while ((entry = readdir(dir)) != NULL) {
        char *path = calloc(strlen(name) + strlen(entry->d_name) + 2, sizeof(char));
        sprintf(path, "%s/%s", name, entry->d_name);
        if ((entry->d_type == DT_DIR) && (recurse != 0)) {
            if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
                free(path);
                continue;
            }
            int rtn = process_dir(path, writer, iter_opts, failures, args, recurse - 1);
            status = max(status, rtn);
        } else {
            // skip
        }
        free(path);
    }
    closedir(dir);

    return status;
}


int process_file(char* fname, fc_writer* writer, const fc_fastq_iter_opts* iter_opts,
        uint64_t* failures, arguments_t* args, int recurse) {
    int status = 0;
    struct stat finfo;
    int res = stat(fname, &finfo);
    if (res == -1) {
        fprintf(stderr, "ERROR  : could not process file %s: %s\n", fname, strerror(errno));
        return errno;
    }

    if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
        if (recurse != 0) {
            char* sfname = strip_trailing_slash(fname);
            int rtn = process_dir(sfname, writer, iter_opts, failures, args, recurse - 1);
            status = max(status, rtn);
            free(sfname);
        }
        return status;
    }

    fc_fastq_iter iter = {0};
    status = fc_fastq_iter_open(&iter, fname, iter_opts);
    if (status != 0) {
        return status;
    }

    fc_writer_begin_file(writer, fname);
    fc_read *fr = NULL;
    while ((fr = fc_fastq_iter_next(&iter)) != NULL) {
        fc_writer_add_read(writer, fr, fname, NULL);
        fc_read_destroy(fr);
    }

    fc_file_stats stats = {0};
    status = fc_fastq_iter_finish(&iter, &stats);
    fc_writer_end_file(writer, fname, &stats);
    if (failures) {
        for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
            failures[i] += stats.failures[i];
        }
    }
    return status;
}


int main(int argc, char **argv) {
    arguments_t args = parse_arguments(argc, argv);

    fc_writer *writer = fc_writer_create();
    if (writer == NULL) exit(1);
    bool use_barcodes = args.demultiplex_dir != NULL;
    size_t barcode_count = use_barcodes ? MAX_BARCODES : 1;
    fc_writer_set_barcode_count(writer, barcode_count, use_barcodes);

    fc_output_sink_opts output_opts = {
        .output_dir = args.demultiplex_dir,
        .reads_per_file = args.reads_per_file,
        .reheader = args.reheader,
        .write_bam = args.write_bam,
        .threads = args.threads
    };
    fc_writer_add_sink(writer, fc_output_sink_create(&output_opts));
    if (args.histograms) {
        fc_stats_engine_config stats_cfg = {
            .track_alignment = false,
            .track_unmapped = true,
            .track_polya = false,
            .polya_cover = 0.0f,
            .polya_qual = 0.0f,
            .polya_rev = false,
            .create_default = !use_barcodes
        };
        fc_writer_configure_stats(writer, &stats_cfg, use_barcodes);

        size_t engine_count = 0;
        fc_stats_engine **stats_engines = fc_writer_get_stats_engines(writer, &engine_count);
        fc_hist_sink_opts hist_opts = {
            .hist_dir = args.histograms,
            .demux_dir = args.demultiplex_dir,
            .use_barcodes = use_barcodes,
            .create_hist_dir = true,
            .write_unmapped = false,
            .write_alignment_metrics = false,
            .write_polya = false,
            .unmapped_is_primary = true,
            .engines = stats_engines,
            .engine_count = engine_count
        };
        fc_writer_add_sink(writer, fc_hist_sink_create(&hist_opts));
    }
    if (args.perread) {
        fc_perread_sink_opts perread_opts = {
            .path = args.perread,
            .sample = args.sample,
            .include_alignment = false,
            .output_unmapped = true
        };
        fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
    }
    if (args.perfile) {
        fc_writer_add_sink(writer, fc_perfile_sink_create(args.perfile, args.sample));
    }
    if (args.runids || args.basecallers) {
        fc_writer_add_sink(writer, fc_counter_sink_create(args.runids, args.basecallers, args.sample));
    }

    fc_fastq_iter_opts iter_opts = {
        .min_length = args.min_length,
        .max_length = args.max_length,
        .min_qscore = args.min_qscore,
        .dust = args.dust,
        .dust_t = args.dust_t,
        .dust_w = args.dust_w,
        .max_dust = args.max_dust,
        .read_group = NULL,
        .tag_name = { 0, 0 },
        .tag_value = 0
    };
    uint64_t failures[NUM_FAILURE_CODES] = {0};

    size_t nfile = 0;
    int status = 0;
    for( ; args.files[nfile] ; nfile++);

    if (nfile==1 && strcmp(args.files[0], "-") == 0) {
        char *ln = NULL;
        size_t n = 0;
        ssize_t nchr = 0;
        int recurse = 0;
        while ((nchr = getline (&ln, &n, stdin)) != -1) {
            ln[strcspn(ln, "\r\n")] = 0;
            int rtn = process_file(ln, writer, &iter_opts, failures, &args, recurse);
            status = max(status, rtn);
        }
        free(ln);
    } else {
        for (size_t i=0; i<nfile; ++i) {
            int rtn = process_file(args.files[i], writer, &iter_opts, failures, &args, args.recurse);
            status = max(status, rtn);
        }
    }

    fc_writer_finish(writer);

    uint64_t total_records =
        failures[R_RECORD_OK]
        + failures[R_TOO_LONG]
        + failures[R_TOO_SHORT]
        + failures[R_LOW_QUALITY]
        + failures[R_DUST_MASKED];
    fprintf(stderr, "INFO   : Processed %" PRIu64 " records in %zu files.\n", total_records, nfile);
    if (status != EXIT_SUCCESS) {
        fprintf(stderr, "WARNING: Error processing files.\n");
    }
    if (!args.force_error) {
        status = EXIT_SUCCESS;
    }

    fprintf(stderr, "\nParsing/filtering summary:\n");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        fprintf(stderr, "%s\t%" PRIu64 "\n", failure_type[i], failures[i]);
    }
    fc_writer_destroy(writer);
    return status;
}
