// bamindex dump program

#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

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

    FILE* fp = fopen(filename, "rb");
    size_t FILE_VERSION = 1;
    size_t written = 0;
    size_t every = 0;
    fread(&FILE_VERSION, sizeof(FILE_VERSION), 1, fp);
    fread(&written, sizeof(written), 1, fp);
    fread(&every, sizeof(every), 1, fp);

    fprintf(stderr, "Reading %zu records from index.\n", written);
    // read: file offset, length qname, qname
    // 4309843968 37 f958cd9c-3d19-43ee-a3ea-aac01c062701
    size_t max_qname = 256;
    char* qname = calloc(max_qname, sizeof(char));
    for (size_t i=0; i<written; ++i) {
        size_t file_offset;
        size_t l_qname;
        fread(&file_offset, sizeof(file_offset), 1, fp);
        fread(&l_qname, sizeof(l_qname), 1, fp);
        if (l_qname > max_qname) {
            qname = realloc(qname, sizeof(char) * l_qname);
            max_qname = l_qname;
        }
        fread(qname, sizeof(char), l_qname, fp);
        fprintf(stdout, "%zu %zu %s\n", file_offset, l_qname, qname);
    }
    free(qname);

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