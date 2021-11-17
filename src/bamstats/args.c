#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <argp.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "args.h"

const char *argp_program_version = "0.4.0";
const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
 "bamstats -- summarise rears/alignments in one or more BAM files.\
 \vThe program creates a simple TSV file containing statistics for \
 each primary alignment stored within the input BAM files.";
static char args_doc[] = "<reference.fasta> <reads.bam> [<reads.bam> ...]";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:"},
    {"region", 'r', "chr:start-end", 0,
        "Genomic region to process."},
    {"threads", 't', "THREADS", 0,
        "Number of threads for BAM processing."},
    {0, 0, 0, 0,
        "Read filtering options:"},
    {"read_group", 'g', "RG", 0,
        "Only process reads from given read group.", 3},
    {"tag_name", 0x100, "TN", 0,
        "Only process reads with a given tag (see --tag_value).", 3},
    {"tag_value", 0x200, "VAL", 0,
        "Only process reads with a given tag value.", 3},
    {"haplotype", 0x300, "VAL", 0,
        "Only process reads from a given haplotype. Equivalent to --tag_name HP --tag_value VAL.", 3},
    { 0 }
};

bool file_exists(char* filename) {
    struct stat st;
    return (stat(filename, &st) == 0);
}

static int tag_items = 0;
static bool tag_given = false;
static bool hp_given = false;
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 'r':
            arguments->region = arg;
            break;
        case 'g':
            arguments->read_group = arg;
            break;
        case 0x100:
            if (strlen(arg) > 2) {
                argp_error(state, "Tag name should be a two-letter code, received: '%s'.", arg);
            }
            memcpy(arguments->tag_name, arg, 2 *sizeof(char));
            tag_items += 1;
            tag_given = true;
            break;
        case 0x200:
            arguments->tag_value = atoi(arg);
            tag_items += 1;
            tag_given = true;
            break;
        case 0x300:
            memcpy(arguments->tag_name, "HP", 2 * sizeof(char));
            arguments->tag_value = atoi(arg);
            tag_items += 2;
            hp_given = true;
            break;
        case 't':
            arguments->threads = atoi(arg);
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                arguments->ref = arg;
                if (!file_exists(arg)) {
                    argp_error(state, "Cannot access reference input file: '%s'.", arg);
                }
                faidx_t *fai = fai_load(arg);
                if (fai == NULL) {
                    argp_error(state, "Cannot read .fasta(.gz) file: '%s'.", arg);
                }
                fai_destroy(fai);
                break;
            } else {
                arguments->bam = (const char**)(&state->argv[state->next - 1]);
                state->next = state->argc;
                break;
            }
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 2)
                argp_usage (state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    args.bam = NULL;
    args.ref = NULL;
    args.region = NULL;
    args.read_group = NULL;
    args.tag_name[0] = '\0';
    args.tag_value = -1;
    args.threads = 1;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    if (tag_items % 2 > 0) {
        fprintf(stderr, "ERROR: Both or neither of --tag_name and --tag_value must be given.\n");
        exit(EXIT_FAILURE);
    }
    if (tag_given && hp_given) {
        fprintf(stderr, "ERROR: If --haplotype is given neither of --tag_name or --tag_value should be provided.\n");
        exit(EXIT_FAILURE);
    }
    return args;
}
