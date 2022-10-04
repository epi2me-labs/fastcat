#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KSEQ_DECLARED

#include "../fastqcomments.h"
#include "args.h"
#include "writer.h"

static inline size_t max ( size_t a, size_t b ) { return a > b ? a : b; }
static inline size_t min ( size_t a, size_t b ) { return a < b ? a : b; }

const double qprobs[100] = {
    1.00000000e+00, 7.94328235e-01, 6.30957344e-01, 5.01187234e-01,
    3.98107171e-01, 3.16227766e-01, 2.51188643e-01, 1.99526231e-01,
    1.58489319e-01, 1.25892541e-01, 1.00000000e-01, 7.94328235e-02,
    6.30957344e-02, 5.01187234e-02, 3.98107171e-02, 3.16227766e-02,
    2.51188643e-02, 1.99526231e-02, 1.58489319e-02, 1.25892541e-02,
    1.00000000e-02, 7.94328235e-03, 6.30957344e-03, 5.01187234e-03,
    3.98107171e-03, 3.16227766e-03, 2.51188643e-03, 1.99526231e-03,
    1.58489319e-03, 1.25892541e-03, 1.00000000e-03, 7.94328235e-04,
    6.30957344e-04, 5.01187234e-04, 3.98107171e-04, 3.16227766e-04,
    2.51188643e-04, 1.99526231e-04, 1.58489319e-04, 1.25892541e-04,
    1.00000000e-04, 7.94328235e-05, 6.30957344e-05, 5.01187234e-05,
    3.98107171e-05, 3.16227766e-05, 2.51188643e-05, 1.99526231e-05,
    1.58489319e-05, 1.25892541e-05, 1.00000000e-05, 7.94328235e-06,
    6.30957344e-06, 5.01187234e-06, 3.98107171e-06, 3.16227766e-06,
    2.51188643e-06, 1.99526231e-06, 1.58489319e-06, 1.25892541e-06,
    1.00000000e-06, 7.94328235e-07, 6.30957344e-07, 5.01187234e-07,
    3.98107171e-07, 3.16227766e-07, 2.51188643e-07, 1.99526231e-07,
    1.58489319e-07, 1.25892541e-07, 1.00000000e-07, 7.94328235e-08,
    6.30957344e-08, 5.01187234e-08, 3.98107171e-08, 3.16227766e-08,
    2.51188643e-08, 1.99526231e-08, 1.58489319e-08, 1.25892541e-08,
    1.00000000e-08, 7.94328235e-09, 6.30957344e-09, 5.01187234e-09,
    3.98107171e-09, 3.16227766e-09, 2.51188643e-09, 1.99526231e-09,
    1.58489319e-09, 1.25892541e-09, 1.00000000e-09, 7.94328235e-10,
    6.30957344e-10, 5.01187234e-10, 3.98107171e-10, 3.16227766e-10,
    2.51188643e-10, 1.99526231e-10, 1.58489319e-10, 1.25892541e-10};


void kahan_sum(double* sum, double term, double* c) {
    double y = term + *c;
    double t = *sum + y;
    *c = (t - *sum) - y;
    *sum = t;
}


float mean_qual(char* qual, size_t len) {
    double qsum = 0;
    double c = 0;
    for (size_t i=0; i<len; ++i) {
        int q = (int)(qual[i]) - 33;
        kahan_sum(&qsum, qprobs[q], &c);
    }
    qsum /= len;
    return -10 * log10(qsum);
}

const char filetypes[4][9] = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};
size_t nfiletypes = 4;

// defined below -- recursion
int process_file(char* fname, writer writer, arguments_t *args);

int process_dir(const char *name, writer writer, arguments_t *args) {
    int status = 0;
    DIR *dir;
    struct dirent *entry;
    char* search;

    if (!(dir = opendir(name))) {
        fprintf(stderr, "Error: could not process directory %s: %s\n", name, strerror(errno));
        return errno;
    }

    while ((entry = readdir(dir)) != NULL) {
        char *path = calloc(strlen(name) + strlen(entry->d_name) + 2, sizeof(char));
        sprintf(path, "%s/%s", name, entry->d_name);
        if (entry->d_type == DT_DIR) {
            if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
                free(path);
                continue;
            }
            int rtn = process_dir(path, writer, args);
            status = max(status, rtn);
        } else {
            for (size_t i=0; i<nfiletypes; ++i) {
                search = strstr(entry->d_name, filetypes[i]);
                if (search != NULL) {
                    fprintf(stderr, "Processing %s\n", path);
                    int rtn = process_file(path, writer, args);
                    status = max(status, rtn);
                    break;
                }
            }
        }
        free(path);
    }
    closedir(dir);
    return status;
}


int process_file(char* fname, writer writer, arguments_t* args) {
    int status = 0;
    struct stat finfo;
    int res = stat(fname, &finfo);
    if (res == -1) {
        // Failure will bubble up through process_file, process_dir to main
        fprintf(stderr, "Error: could not process file %s: %s\n", fname, strerror(errno));
        return errno;
    }
    if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
        if (args->recurse) {
            char* sfname = strip_path(fname);
            int rtn = process_dir(sfname, writer, args);
            status = max(status, rtn);
            free(sfname);
        } else {
            fprintf(stderr, "Warning: input '%s' is a directory and -x (recursive mode) was not given.\n", fname);
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
    while ((status = kseq_read(seq)) >= 0) {
        ++n ; slen += seq->seq.l;
        minl = min(minl, seq->seq.l);
        maxl = max(maxl, seq->seq.l);
        if (seq->qual.l == 0) { status = -99; break; }
        float mean_q = mean_qual(seq->qual.s, seq->qual.l);
        kahan_sum(&meanq, mean_q, &c);
        read_meta meta = parse_read_meta(seq->comment);
        if ((seq->seq.l >= args->min_length) && (seq->seq.l <= args->max_length) && (mean_q >= args->min_qscore)) {
            write_read(writer, seq, meta, mean_q, fname);
        }
        destroy_read_meta(meta);
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
        fprintf(writer->perfile, "%zu\t%zu\t%zu\t%zu\t%1.2f\n", n, slen, minl, maxl, meanq/n);
    }

    // cleanup
    kseq_destroy(seq);
    gzclose(fp);
    return status;
}


int main(int argc, char **argv) {
    arguments_t args = parse_arguments(argc, argv);
    writer writer = initialize_writer(args.demultiplex_dir, args.perread, args.perfile, args.sample, args.reheader);
    if (writer == NULL) exit(1);

    int nfile = 0;
    int status = 0;
    for( ; args.files[nfile] ; nfile++);

    if (nfile==1 && strcmp(args.files[0], "-") == 0) {
        char *ln = NULL;
        size_t n = 0;
        ssize_t nchr = 0;
        while ((nchr = getline (&ln, &n, stdin)) != -1) {
            ln[strcspn(ln, "\r\n")] = 0;
            int rtn = process_file(ln, writer, &args);
            status = max(status, rtn);
        }
        free(ln);
    } else {
        for (size_t i=0; i<nfile; ++i) {
            int rtn = process_file(args.files[i], writer, &args);
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
