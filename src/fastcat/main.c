#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KSEQ_DECLARED

#include "../common.h"
#include "../fastqcomments.h"
#include "../kh_counter.h"
#include "../sdust/sdust.h"
#include "args.h"
#include "parsing.h"
#include "writer.h"
#include "../fc_read.h"


const char filetypes[4][10] = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};
size_t nfiletypes = 4;

// defined below -- recursion
int process_file(char* fname, writer writer, arguments_t *args, int recurse);

int process_dir(const char *name, writer writer, arguments_t *args, int recurse) {
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
                    int rtn = process_file(path, writer, args, recurse - 1);
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
            int rtn = process_dir(path, writer, args, recurse - 1);
            status = max(status, rtn);
        } else {
            // skip
        }
        free(path);
    }
    closedir(dir);

    return status;
}


double dust_fraction(uint8_t *seq, size_t len, int t, int w) {
    if (len == 0) return 0.0;
    uint64_t *r;
    int n;
    r = sdust(0, (uint8_t*) seq, -1, t, w, &n);
    int masked_bases = 0;
    for (int i = 0; i < n; ++i) {
        int start = (int)(r[i] >> 32);
        int end = (int)r[i];
        masked_bases += (end - start);
    }
    return (double) masked_bases / len;
}


int process_file(char* fname, writer writer, arguments_t* args, int recurse) {
    int status = 0;
    struct stat finfo;
    int res = stat(fname, &finfo);
    if (res == -1) {
        fprintf(stderr, "ERROR  : could not process file %s: %s\n", fname, strerror(errno));
        return errno;
    }

    // handle directory input
    if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
        if (recurse != 0) {
            char* sfname = strip_path(fname);
            int rtn = process_dir(sfname, writer, args, recurse - 1);
            status = max(status, rtn);
            free(sfname);
        }
        return status;
    }

    // these are all per-file stats, separate from writer stats
    size_t n_valid = 0, total_bases = 0;
    size_t min_len = UINTMAX_MAX, max_len = 0;
    double average_q = 0.0, q_sum = 0.0;
    kh_counter_t *run_ids = kh_counter_init();
    kh_counter_t *basecallers = kh_counter_init();
    uint64_t failures[NUM_FAILURE_CODES] = {0};

    gzFile fp;
    kseq_t *seq;
    fp = gzopen(fname, "r");
    seq = kseq_init(fp);
    bool truncated = false;  // track if last read record was truncated
    status = 0;
    while ((status = kseq_read(seq)) != -1) {  // EOF - normal exit
        if (status == -2) {  // truncated quality string
            failures[F_QUAL_TRUNCATED]++;
            truncated = true;
            continue;
        } else if (status == -3) {  // error reading stream
            failures[F_STREAM_ERROR]++;
            break;
        } else if (status < 0) {  // other errors
            failures[F_UNKNOWN_ERROR]++;
            break;
        }
        if (seq->qual.l == 0) {
            failures[F_QUAL_MISSING]++;
            truncated = true; // not present is truncated \:D/
            status = -2;
            continue;
        } else {
            truncated = false;
            failures[R_RECORD_OK]++;
        }

        // accumulate stats only for reads within length and quality thresholds
        if (seq->seq.l > args->max_length) {
            failures[R_TOO_LONG]++;
            continue;
        }
        if (seq->seq.l < args->min_length) {
            failures[R_TOO_SHORT]++;
            continue;
        }

        // build fc_read now as it handles mean_q from meta or recomputing
        fc_read* fr = fc_read_from_kseq(seq);
        float mean_q = fr->meta->mean_q;
        if (mean_q < args->min_qscore) {
            fc_read_destroy(fr);
            failures[R_LOW_QUALITY]++;
            continue;
        }

        if (args->dust) {
            double masked_fraction = dust_fraction((uint8_t*)seq->seq.s, seq->seq.l, args->dust_t, args->dust_w);
            if (masked_fraction > args->max_dust) {
                failures[R_DUST_MASKED]++;
                continue;
            }
        }

        // update per-file stats
        ++n_valid ; total_bases += seq->seq.l;
        min_len = min(min_len, seq->seq.l);
        max_len = max(max_len, seq->seq.l);
        kahan_sum(&average_q, mean_q, &q_sum);
        kh_counter_increment(run_ids, fr->meta->runid);
        kh_counter_increment(basecallers, fr->meta->basecaller);

        // write read to output
        write_read(writer, fr, fname);
        fc_read_destroy(fr);
    }

    // handle final reader status
    if (truncated) {
        // truncated treated as stream error
        failures[F_STREAM_ERROR]++;
        status = -3;
    }
    else if (status == -1) {
        failures[F_FILE_OK]++;
    }

    status = status == -1 ? EXIT_SUCCESS : EXIT_FAILURE;

    if (failures[F_STREAM_ERROR] > 0) {
        fprintf(stderr,
            "WARNING: file '%s' is possibly truncated.\n",
            fname);
    } else if (failures[F_QUAL_MISSING] > 0) {
        fprintf(stderr,
           "WARNING: no quality string found for %" PRIu64 " records in file '%s'.\n",
           failures[F_QUAL_MISSING], fname);
    } else if (failures[F_QUAL_TRUNCATED] > 0) {
        fprintf(stderr,
            "WARNING: truncated quality string found for %" PRIu64 " records in file '%s'.\n",
            failures[F_QUAL_TRUNCATED], fname);
    } else if (failures[F_UNKNOWN_ERROR] > 0) {
        fprintf(stderr, "WARNING: unknown error reading file '%s'.\n", fname);
    }

    // summary entries
    // TODO: move this to writer, we can wrap per file stats in a struct and pass that
    if(writer->perfile != NULL) {
        fprintf(writer->perfile, "%s\t", fname);
        if (writer->sample != NULL) fprintf(writer->perfile, "%s\t", args->sample);
        if (n_valid == 0) {
            // there were no reads in the input file
            fprintf(writer->perfile, "0\t0\t0\t0\t0.00");
        } else {
            fprintf(writer->perfile, "%zu\t%zu\t%zu\t%zu\t%.2f",
                n_valid, total_bases, min_len, max_len, average_q / n_valid
            );
        }
        for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
            fprintf(writer->perfile, "\t%" PRIu64, failures[i]);
        }
        fprintf(writer->perfile, "\n");
    }
    if(writer->runids != NULL) {
        for (khiter_t k = 0; k < kh_end(run_ids); ++k) {
            if (kh_exist(run_ids, k)) {
                fprintf(writer->runids, "%s\t", fname);
                if (writer->sample != NULL) fprintf(writer->runids, "%s\t", args->sample);
                fprintf(writer->runids, "%s\t%d\n", kh_key(run_ids, k), kh_val(run_ids, k));
            }
        }
    }
    if(writer->basecallers != NULL) {
        for (khiter_t k = 0; k < kh_end(basecallers); ++k) {
            if (kh_exist(basecallers, k)) {
                fprintf(writer->basecallers, "%s\t", fname);
                if (writer->sample != NULL) fprintf(writer->basecallers, "%s\t", args->sample);
                fprintf(writer->basecallers, "%s\t%d\n", kh_key(basecallers, k), kh_val(basecallers, k));
            }
        }
    }
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        writer->failures[i] += failures[i];
    }

    // cleanup
    kh_counter_destroy(basecallers);
    kh_counter_destroy(run_ids);
    kseq_destroy(seq);
    gzclose(fp);
    return status;
}


int main(int argc, char **argv) {
    arguments_t args = parse_arguments(argc, argv);

    writer writer = initialize_writer(
        args.demultiplex_dir, args.histograms, args.perread, args.perfile,
        args.runids, args.basecallers, args.sample,
        args.reheader, args.write_bam, args.reads_per_file,
        args.threads);
    if (writer == NULL) exit(1);

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
            int rtn = process_file(ln, writer, &args, recurse);
            status = max(status, rtn);
        }
        free(ln);
    } else {
        for (size_t i=0; i<nfile; ++i) {
            int rtn = process_file(args.files[i], writer, &args, args.recurse);
            status = max(status, rtn);
        }
    }

    uint64_t total_records =
        writer->failures[R_RECORD_OK] 
        + writer->failures[R_TOO_LONG]
        + writer->failures[R_TOO_SHORT]
        + writer->failures[R_LOW_QUALITY]
        + writer->failures[R_DUST_MASKED];
    fprintf(stderr, "INFO   : Processed %" PRIu64 " records in %zu files.\n", total_records, nfile);
    if (status != EXIT_SUCCESS) {
        fprintf(stderr, "WARNING: Error processing files.\n");
    }
    if (!args.force_error) {
        status = EXIT_SUCCESS;
    }

    fprintf(stderr, "\nParsing/filtering summary:\n");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        fprintf(stderr, "%s\t%" PRIu64 "\n", failure_type[i], writer->failures[i]);
    }
    destroy_writer(writer);
    return status;
}
