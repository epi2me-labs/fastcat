// bamindex fetch program

#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "htslib/bgzf.h"

#include "common.h"

#include <argp.h>

static char doc[] = 
"bamindex fetch -- fetch records from a BAM according to an index.\
\vThe program simply will fetch a batch of records from a BAM file" \
"using and index and a chunk ID";
static char args_doc[] = "<reads.bam.bci>";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:"},
    {"threads", 't', "THREADS", 0,
        "Number of threads for BAM processing."},
    {"chunk", 'c', "SIZE", 0,
        "Chunk index to retrieve."},
    { 0 }
};

typedef struct arguments {
    const char* bam;
    const char* index;
    int chunk_idx;
    int threads;
} arguments_t;

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 't':
            arguments->threads = atoi(arg);
            break;
        case 'c':
            arguments->chunk_idx = atoi(arg);
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
    args.index = NULL;
    args.chunk_idx = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}


void index_fetch(const char* bam_fname, const char* index_fname, int chunk, int threads) {
    htsFile *fp = hts_open(bam_fname, "r");
    bam_hdr_t *h = sam_hdr_read(fp);
    if(fp == NULL || h == NULL) {
        fprintf(stderr, "Could not open %s\n", bam_fname);
        exit(EXIT_FAILURE);
    }
    htsThreadPool p = {NULL, 0};
    if (threads > 1 ) {
        p.pool = hts_tpool_init(threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    struct stat st;
    if (stat(index_fname, &st) != 0) {
        errx(1, "Cannot open index file %s\n", index_fname);
        exit(EXIT_FAILURE);
    }
    FILE* idx_fp = fopen(index_fname, "rb");
    size_t FILE_VERSION = 1;
    size_t written = 0;
    size_t chunk_size = 0;
    fread(&FILE_VERSION, sizeof(FILE_VERSION), 1, idx_fp);
    fread(&written, sizeof(written), 1, idx_fp);
    fread(&chunk_size, sizeof(chunk_size), 1, idx_fp);

    // read: file offset, length qname, qname
    // 4309843968 37 f958cd9c-3d19-43ee-a3ea-aac01c062701
    size_t max_qname = 256;
    char* qname = calloc(max_qname, sizeof(char));
    size_t file_offset;
    size_t l_qname;
    int i = 0;
    while (i <= chunk) {
        fread(&file_offset, sizeof(file_offset), 1, idx_fp);
        fread(&l_qname, sizeof(l_qname), 1, idx_fp);
        if (l_qname > max_qname) {
            qname = realloc(qname, sizeof(char) * l_qname);
            max_qname = l_qname;
        }
        fread(qname, sizeof(char), l_qname, idx_fp);
        i++;
    }
    fprintf(stderr, "Starting from: %zu %zu %s\n", file_offset, l_qname, qname);
    fprintf(stderr, "Reading %zu records from bam.\n", chunk_size);


    if(bgzf_seek(fp->fp.bgzf, file_offset, SEEK_SET) != 0) {
        fprintf(stderr, "Failed to seek to first record.\n");
        exit(EXIT_FAILURE);
    }

    int ret = 0;
    written = 0;
    bam1_t* b = bam_init1();
    htsFile * out_fp;
    if ((out_fp = hts_open("-", "wb0")) == 0) {
        fprintf(stderr, "Failed to open standard output for writing.\n");
        exit(EXIT_FAILURE);
    }

    // TODO: fill in the NULLs here
    if(sam_hdr_add_pg(h, "bamindex.fetch", "VN", argp_program_version, NULL, NULL, NULL) != 0){
        fprintf(stderr, "Failed to add PG line to the header.\n");
        exit(EXIT_FAILURE);
    }
    if(sam_hdr_write(out_fp, h) != 0) {
        fprintf(stderr, "Failed to write the SAM header.\n");
        exit(EXIT_FAILURE);
    }
    while ((sam_read1(fp, h, b) >= 0) && (written < (int)chunk_size)) {
        if((sam_write1(out_fp, h, b) < 0)) {
            fprintf(stderr, "Failed to write output record.");
            exit(EXIT_FAILURE);
        }
        written++;
    }
    hts_close(out_fp);

    free(qname);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    fprintf(stderr, "Written %zu records to output.\n", written);
}


int main_fetch(int argc, char *argv[]) {
    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
    args.index = generate_index_filename(args.bam, args.index);
    index_fetch(args.bam, args.index, args.chunk_idx, args.threads);
    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    return EXIT_SUCCESS;
}