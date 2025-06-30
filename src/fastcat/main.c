#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KSEQ_DECLARED

#include "../common.h"
#include "../fastqcomments.h"
#include "../kh_counter.h"
#include "../sdust/sdust.h"
#include "args.h"
#include "writer.h"


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
        fprintf(stderr, "Error: could not process directory %s: %s\n", name, strerror(errno));
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
        fprintf(stderr, "Error: could not process directory %s: %s\n", name, strerror(errno));
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
        fprintf(stderr, "Error: could not process file %s: %s\n", fname, strerror(errno));
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

    gzFile fp;
    kseq_t *seq;

    fp = gzopen(fname, "r");
    seq = kseq_init(fp);
    size_t n = 0, slen = 0;
    size_t minl = UINTMAX_MAX, maxl = 0;
    double meanq = 0.0, c = 0.0;
    status = 0;
    kh_counter_t *run_ids = kh_counter_init();
    kh_counter_t *basecallers = kh_counter_init();
    while ((status = kseq_read(seq)) >= 0) {
        // accumulate stats only for reads within length and quality thresholds
        if (seq->qual.l == 0) { status = -99; break; }
        if ((seq->seq.l >= args->min_length) && (seq->seq.l <= args->max_length)) {
            float mean_q = mean_qual_naive(seq->qual.s, seq->qual.l);
            if (mean_q < args->min_qscore) continue;

            // do some housework
            if (args->dust) {
                double masked_fraction = dust_fraction((uint8_t*)seq->seq.s, seq->seq.l, args->dust_t, args->dust_w);
                if (masked_fraction > args->max_dust) {
                    if (args->verbose) {
                        fprintf(stderr, "Skipping read %s due to excessive dust masking (%.2f > %.2f)\n",
                            seq->name.s, masked_fraction, args->max_dust);
                    }
                    continue;
                }
            }

            ++n ; slen += seq->seq.l;
            minl = min(minl, seq->seq.l);
            maxl = max(maxl, seq->seq.l);
            kahan_sum(&meanq, mean_q, &c);
            read_meta meta = parse_read_meta(seq->comment);
            write_read(writer, seq, meta, mean_q, fname);
            kh_counter_increment(run_ids, meta->runid);
            kh_counter_increment(basecallers, meta->basecaller);
            destroy_read_meta(meta);
        }
    }

    // handle errors
    switch (status) {
        case -1:
            status = EXIT_SUCCESS;
            break;
        case -2:
            status = EXIT_FAILURE;
            fprintf(stderr, "Truncated quality string found for record in file '%s'.\n", fname);
            break;
        case -3:
            status = EXIT_FAILURE;
            fprintf(stderr, "Error reading file '%s', possibly truncated\n", fname);
            break;
        case -99:
            status = EXIT_FAILURE;
            fprintf(stderr, "No quality string found for record in file '%s' (FASTA is unsupported).\n", fname);
            break;
        default:
            status = EXIT_FAILURE;
            fprintf(stderr, "Unknown error reading file '%s'.\n", fname);
    }

    // summary entries
    if(writer->perfile != NULL) {
        fprintf(writer->perfile, "%s\t", fname);
        if (writer->sample != NULL) fprintf(writer->perfile, "%s\t", args->sample);
        if (n == 0) {
            // there were no reads in the input file
            fprintf(writer->perfile, "0\t0\t0\t0\t0.00\n");
        } else {
            fprintf(writer->perfile, "%zu\t%zu\t%zu\t%zu\t%.2f\n",
                n, slen, minl, maxl, meanq/n
            );
        }
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
    destroy_writer(writer);

    if (status != 0) {
        fprintf(stderr, "Completed processing with errors. Outputs may be incomplete.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
