// bamindex build program

#include <err.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "htslib/bgzf.h"

#include "common.h"
#include "../version.h"

#include <argp.h>

typedef struct arguments {
    const char* bam;
    int threads;
    int chunk_size;
} arguments_t;

static char doc[] = 
"bamindex build -- create a BAM index corresponding to batches of records.\
\vThe program creates a simple index of file offsets for each of every \
(n * M)th alignment record. No care is taken to keep records corresponding \
to the same query together, or any other such niceities. Its intended to \
be used simply with unaligned, unsorted BAMs.";
static char args_doc[] = "<reads.bam>";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:"},
    {"threads", 't', "THREADS", 0,
        "Number of threads for BAM processing."},
    {"chunk_size", 'c', "SIZE", 0,
        "Number of records in a chunk."},
    { 0 }
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 't':
            arguments->threads = atoi(arg);
            break;
        case 'c':
            arguments->chunk_size = atoi(arg);
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                arguments->bam = arg;
                break;
            }
            break;
        case ARGP_KEY_END:
            if (state->arg_num != 1)
                argp_usage (state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

static arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    args.bam = NULL;
    args.threads = 1;
    args.chunk_size = 1;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}

void index_build(const char* filename, const char* output_fname, int threads, size_t every) {
    htsFile *fp = hts_open(filename, "r");
    bam_hdr_t *h = sam_hdr_read(fp);
    if(fp == NULL || h == NULL) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(EXIT_FAILURE);
    }
    
    FILE* out_fp = fopen(output_fname, "wb");
    size_t FILE_VERSION = 1;
    size_t written = 0;
    fwrite(&FILE_VERSION, sizeof(FILE_VERSION), 1, out_fp);
    fwrite(&written, sizeof(written), 1, out_fp);
    fwrite(&every, sizeof(every), 1, out_fp);

    htsThreadPool p = {NULL, 0};
    if (threads > 1 ) {
        p.pool = hts_tpool_init(threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    int ret = 0;
    int idx = 0;
    bam1_t* b = bam_init1();
    size_t file_offset = bgzf_tell(fp->fp.bgzf);
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        if ((idx % every) != 0) {
            file_offset = bgzf_tell(fp->fp.bgzf);
            idx++;
            continue;
        }
        char* qname = bam_get_qname(b);
        size_t l_qname = strlen(qname) + 1;
        written++;
        // write: file offset, length qname, qname
        fwrite(&file_offset, sizeof(file_offset), 1, out_fp);
        fwrite(&l_qname, sizeof(l_qname), 1, out_fp);
        fwrite(qname, sizeof(char), l_qname, out_fp);

        if(idx % 100000 == 0) {
            fprintf(stderr, "Record %d %zu\n", idx, file_offset);
        }
        file_offset = bgzf_tell(fp->fp.bgzf);
        idx++;
    }

    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    // fill in how many records we wrote
    fseek(out_fp, sizeof(FILE_VERSION), SEEK_SET);
    fwrite(&written, sizeof(written), 1, out_fp);
    fclose(out_fp);
    fprintf(stderr, "Written %zu/%d records to index.\n", written, idx+1);
}


int main_build(int argc, char *argv[]) {
    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
#ifdef NOTHREADS
    if (args.threads != 1) {
        fprintf(
            stderr,
            "--threads set to %d, but threading not supported by this build.\n", args.threads);
    }
#endif

    char* index_fname = generate_index_filename(args.bam, NULL);
    fprintf(stderr, "threads %d\n", args.threads);
    fprintf(stderr, "chunk_size %d\n", args.chunk_size);
    index_build(args.bam, index_fname, args.threads, args.chunk_size);
    free(index_fname);

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    return EXIT_SUCCESS;
}
