#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <argp.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "args.h"
#include "common.h"
#include "../version.h"

const char *argp_program_bug_address = "support@nanoporetech.com";
static char doc[] = 
"bamstats -- summarise reads/alignments in and input BAM file.\
\vThe program creates a simple TSV file containing statistics for \
each primary alignment stored within the input BAM file.";
static char args_doc[] = "<reads.bam>";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:", 0},
    {"region", 'r', "chr:start-end", 0,
        "Genomic region to process.", 0},
    {"bed", 'b', "BEDFILE", 0,
        "BED file for regions to process.", 0},
    {"threads", 't', "THREADS", 0,
        "Number of threads (for BAM decompression and BED output compression).", 0},
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
        "Include unmapped/unplaced reads in output.", 2},
    {"read_group", 'g', "RG", 0,
        "Only process reads from given read group.", 2},
    {"tag_name", 0x100, "TN", 0,
        "Only process reads with a given tag (see --tag_value).", 2},
    {"tag_value", 0x200, "VAL", 0,
        "Only process reads with a given tag value.", 2},
    {"haplotype", 0x300, "VAL", 0,
        "Only process reads from a given haplotype. Equivalent to --tag_name HP --tag_value VAL.", 2},

    {0, 0, 0, 0,
        "Coverage calculation options:", 0},
    {0, 0, 0, 0,
        "Outputs produced are similar to that by mosdepth. Outputs for all sequences in the BAM index are produced in the top-level directory (stratified by reference sequence). An identical set of files is produced (in sub-directories) for each BED file provided. Each BED file must have a corresponding name provided via --coverage_names, which is used to name the sub-directory. The --segments option is used to produce outputs across the full reference broken into fixed-length segments, and output to segments-LENGTH. The use of the segments option with small lengths can significantly affect the performance of the program: it should be preferred to use the per-base output and segment after the fact.\n" , 3},
    {"coverage", 0x1004, "DIRECTORY", 0,
        "Enable coverage calculations and output to provided directory. (default: bamstats-coverages)", 3},
    {"coverage_beds", 0x1000, "BEDFILE ...", 0,
        "BED file(s) for calculating coverage (space-separated list).", 3},
    {"coverage_names", 0x1001, "NAME ...", 0,
        "Name(s) for coverage BED file(s) (space-separated list).", 3},
    {"thresholds", 0x1002, "VALUES ...", 0,
        "Coverage thresholds to produce sparse cumulative distribution (space-separated integers).", 3},
    {"segments", 0x1003, "LENGTH ...", 0,
        "Segment length(s) for which to produce outputs. (space-separated integers).", 3},

    {0, 0, 0, 0,
        "Poly-A Options:", 0},
    {"poly_a", 0x500, 0, 0,
        "Enable poly-A tail length histogram.", 4},
    {"poly_a_cover", 0x600, "PCT_COVERAGE", 0,
        "Reference alignment coverage for acceptance of read. (default: 95)", 4},
    {"poly_a_qual", 0x700, "QUAL", 0,
        "Read mean Q score for acceptance of read. (default: 10)", 4},
    {"poly_a_rev", 0x800, 0, 0,
        "Allow reverse alignments (useful for cDNA, default is appropriate for direct RNA seq).", 4},
    { 0 }
};


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
        case 0x1000:
            slurp_args(&arguments->coverage_beds, &arguments->n_coverage_beds, arg, state);
            break;
        case 0x1001:
            slurp_args(&arguments->coverage_names, &arguments->n_coverage_names, arg, state);
            break;
        case 0x1002:
            slurp_ints(&arguments->coverage_thresholds, &arguments->n_coverage_thresholds, arg, state);
            break;
        case 0x1003:
            slurp_ints(&arguments->segments, &arguments->n_segments, arg, state);
            break;
        case 0x1004:
            arguments->coverage = true;
            arguments->coverages = arg;
            break;

        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                arguments->bam = arg;
                break;
            }
            argp_usage(state);
            break;

        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;

        case ARGP_KEY_END:
            if (!arguments->bam) argp_error(state, "Missing <reads.bam>");
            if (arguments->n_coverage_beds && arguments->n_coverage_beds != arguments->n_coverage_names) {
                argp_error(state, "Mismatched counts: --coverage_beds (%zu) vs --coverage_names (%zu).",
                           arguments->n_coverage_beds, arguments->n_coverage_names);
            }
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
    args.coverage = false;
    args.coverages = "bamstats-coverages";  // will be overwrite by user
    args.coverage_beds = NULL;
    args.n_coverage_beds = 0;
    args.coverage_names = NULL;
    args.n_coverage_names = 0;
    args.coverage_thresholds = NULL;
    args.n_coverage_thresholds = 0;
    args.segments = NULL;
    args.n_segments = 0;
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

void destroy_args(arguments_t *args) {
    if (args->coverage_beds) {
        free(args->coverage_beds);
    }
    if (args->coverage_names) {
        free(args->coverage_names);
    }
    if (args->coverage_thresholds) {
        free(args->coverage_thresholds);
    }
    if (args->segments) {
        free(args->segments);
    }
}
