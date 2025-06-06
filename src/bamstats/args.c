#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <argp.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "args.h"
#include "../version.h"

const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
"bamstats -- summarise rears/alignments in one or more BAM files.\
\vThe program creates a simple TSV file containing statistics for \
each primary alignment stored within the input BAM files.";
static char args_doc[] = "<reads.bam>";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:", 0},
    {"region", 'r', "chr:start-end", 0,
        "Genomic region to process.", 0},
    {"bed", 'b', "BEDFILE", 0,
        "BED file for regions to process.", 0},
    {"threads", 't', "THREADS", 0,
        "Number of threads for BAM processing.", 0},
    {"sample", 's',"SAMPLE NAME",   0,
        "Sample name (if given, adds a 'sample_name' column).", 0},
    {"flagstats", 'f', "FLAGSTATS", 0,
        "File for outputting alignment flag counts.", 0},
    {"runids", 'i', "ID SUMMARY",  0,
        "Run ID summary output", 0},
    {"basecallers", 'l', "BASECALLERS", 0,
        "Basecaller summary output", 0},
    {"histograms", 0x400, "DIRECTORY", 0,
        "Directory for outputting histogram information. (default: bamstats-histograms)", 0},
    {"recalc_qual", 0x900, 0, 0,
        "Force recomputing mean quality, else use 'qs' tag in BAM if present.", 0},
    {0, 0, 0, 0,
        "Read filtering options:", 0},
    {"unmapped", 'u', 0, 0,
        "Include unmapped/unplaced reads in output.", 3},
    {"read_group", 'g', "RG", 0,
        "Only process reads from given read group.", 3},
    {"tag_name", 0x100, "TN", 0,
        "Only process reads with a given tag (see --tag_value).", 3},
    {"tag_value", 0x200, "VAL", 0,
        "Only process reads with a given tag value.", 3},
    {"haplotype", 0x300, "VAL", 0,
        "Only process reads from a given haplotype. Equivalent to --tag_name HP --tag_value VAL.", 3},
    {0, 0, 0, 0,
        "Poly-A Options:", 0},
    {"poly_a", 0x500, 0, 0,
        "Enable poly-A tail length histogram.", 5},
    {"poly_a_cover", 0x600, "PCT_COVERAGE", 0,
        "Reference alignment coverage for acceptance of read. (default: 95)", 5},
    {"poly_a_qual", 0x700, "QUAL", 0,
        "Read mean Q score for acceptance of read. (default: 10)", 5},
    {"poly_a_rev", 0x800, 0, 0,
        "Allow reverse alignments (useful for cDNA, default is appropriate for direct RNA seq).", 5},
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
        case 'b':
            arguments->bed = arg;
            break;
        case 'g':
            arguments->read_group = arg;
            break;
        case 'f':
            arguments->flagstats = arg;
            break;
        case 'i':
            arguments->runids = arg;
            break;
        case 'l':
            arguments->basecallers = arg;
            break;
        case 0x400:
            arguments->histograms = arg;
            break;
        case 0x500:
            arguments->poly_a = true;
            break;
        case 0x600:
            arguments->poly_a_cover = atof(arg);
            break;
        case 0x700:
            arguments->poly_a_qual = atof(arg);
            break;
        case 0x800:
            arguments->poly_a_rev = true;
            break;
        case 's':
            arguments->sample = arg;
            break;
        case 'u':
            arguments->unmapped = true;
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
        case 0x900:
            arguments->force_recalc_qual = true;
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                arguments->bam = (const char**)(&state->argv[state->next - 1]);
                state->next = state->argc;
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

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    args.bam = NULL;
    args.flagstats = NULL;
    args.runids = NULL;
    args.basecallers = NULL;
    args.histograms = "bamstats-histograms";
    args.poly_a = false;
    args.poly_a_cover = 95;
    args.poly_a_qual = 10;
    args.poly_a_rev = false;
    args.sample = NULL;
    args.ref = NULL;
    args.region = NULL;
    args.bed = NULL;
    args.unmapped = false;
    args.read_group = NULL;
    args.tag_name[0] = '\0';
    args.tag_value = -1;
    args.threads = 1;
    args.force_recalc_qual = false;
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
