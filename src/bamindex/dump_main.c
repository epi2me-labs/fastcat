// bamindex dump program

#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "index.h"

#include <argp.h>

static char doc[] = 
"bamindex dump -- dump a BAM chunk index to stdout as text.\
\vThe program simply writes the contents of an index to stdout for human \
inspection. It has no other purpose.";
static char args_doc[] = "<reads.bam.bci>";
static struct argp_option options[] = {
    { 0 }
};

typedef struct arguments {
    const char* index;
} arguments_t;

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                arguments->index = arg;
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
    args.index = NULL;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}


void index_dump(const char* filename) {
    struct stat st;
    if (stat(filename, &st) != 0) {
        errx(1, "Cannot open index file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(filename, "rb");
    bc_idx_t *idx;
    if((idx = bc_idx_read(fp)) == NULL) {
        fprintf(stderr, "Couldn't read index file: %s.\n", filename);
        exit(EXIT_FAILURE);
    }
    fprintf(stdout, "Index contains %zu chunks of size %zu.\n", idx->n_chunks, idx->chunk_size);
    for (size_t i=0; i<idx->n_chunks; ++i){
        fprintf(stdout, "%zu %s\n", (idx->recs[i]).file_offset, (idx->recs[i]).qname);
    }
    bc_idx_destroy(idx);
    fclose(fp);
}


int main_dump(int argc, char *argv[]) {
    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
    index_dump(args.index);
    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    return EXIT_SUCCESS;
}