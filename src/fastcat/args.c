#include <stdlib.h>
#include <argp.h>

#include "args.h"
#include "../version.h"

const char *argp_program_bug_address = "support@nanoporetech.com";
static char doc[] = 
"fastcat -- concatenate and summarise .fastq(.gz) files.\
\vInput files may be given on stdin by specifing the input as '-'. \
Also accepts directories as input and looks for .fastq(.gz) files in \
the top-level directory. Recurses into sub-directories when the \
-x option is given. The command \
will exit non-zero if any file encountered cannot be read.";
static char args_doc[] = "reads1.fastq(.gz) reads2.fastq(.gz) dir-with-fastq ...";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:", 0},
    {"recurse", 'x', 0, 0,
        "Search directories recursively for '.fastq', '.fq', '.fastq.gz', and '.fq.gz' files.", 0},
    {"threads", 't', "THREADS", 0,
        "Number of threads for output compression (only with --bam_out.", 0},
    {"force_error", 'e', 0, 0,
        "Exit with non-zero status if any files, or records, contained errors.", 0},

    {0, 0, 0, 0,
        "Output options:", 0},
    {"sample", 's', "SAMPLE NAME",   0,
        "Sample name (if given, adds a 'sample_name' column).", 0},
    {"reads_per_file", 'c', "NUM", 0,
        "Split reads into files with a set number of reads (default: single file).", 0},
    {"reheader", 'H', 0, 0,
        "Rewrite fastq header comments as SAM tags (useful for passing through minimap2).", 0},
    {"bam_out", 'B', 0, 0,
        "Output data as unaligned BAM.", 0},
    {"verbose", 'v', 0, 0,
        "Verbose output.", 0},

    {0, 0, 0, 0,
        "Output file selection:", 0},
    {"read", 'r', "READ SUMMARY",  0,
        "Per-read summary output", 0},
    {"file", 'f', "FILE SUMMARY",  0,
        "Per-file summary output", 0},
    {"runids", 'i', "ID SUMMARY",  0,
        "Run ID summary output", 0},
    {"basecallers", 'l', "CALLER SUMMARY",  0,
        "Basecaller mode summary output", 0},
    {"demultiplex", 'd', "OUT DIR",  0,
        "Separate barcoded samples using fastq header information. Option value is top-level output directory.", 0},
    {"histograms", 0x400, "DIRECTORY", 0,
        "Directory for outputting histogram information. When --demultiplex is enabled histograms are written to per-sample demultiplexed output directories. (default: fastcat-histograms)", 0},

    {0, 0, 0, 0,
        "Read filtering options:", 0},
    {"min_length", 'a', "MIN READ LENGTH", 0,
        "minimum read length to output (excluded reads remain listed in summaries).", 0},
    {"max_length", 'b', "MAX READ LENGTH", 0,
        "maximum read length to output (excluded reads remain listed in summaries).", 0},
    {"min_qscore", 'q', "MIN READ QSCORE", 0,
        "minimum read Qscore to output (excluded reads remain listed in summaries).", 0},
    {"dust", 0x500, 0, 0,
        "Enable DUST filtering of reads (default: disabled).", 0},

    { 0, 0, 0, 0,
        "Advanced cleaning options:", 0},
    {"max_dust", 0x600, "MAX DUST", 0,
        "Maximum proportion of low-complexity regions to allow in reads (default: 0.95).", 0},
    {"dust_w", 0x700, "DUST W", 0,
        "Window size for DUST filtering (default: 64).", 0},
    {"dust_t", 0x800, "DUST T", 0,
        "Threshold for DUST filtering (default: 20).", 0},
    { 0 }
};


static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 'r':
            arguments->perread = arg;
            break;
        case 'f':
            arguments->perfile = arg;
            break;
        case 'i':
            arguments->runids = arg;
            break;
        case 'l':
            arguments->basecallers = arg;
            break;
        case 's':
            arguments->sample = arg;
            break;
        case 'a':
            arguments->min_length = atoi(arg);
            break;
        case 'b':
            arguments->max_length = atoi(arg);
            break;
        case 'c':
            arguments->reads_per_file = atoi(arg);
            break;
        case 'd':
            arguments->demultiplex_dir = arg;
            break;
        case 0x400:
            arguments->histograms = arg;
            break;
        case 0x500:
            arguments->dust = 1;
            break;
        case 0x600:
            arguments->max_dust = atof(arg);
            if (arguments->max_dust < 0 || arguments->max_dust > 1) {
                argp_error(state, "max_dust must be between 0 and 1.");
            }
            break;
        case 0x700:
            arguments->dust_w = atoi(arg);
            if (arguments->dust_w <= 0) {
                argp_error(state, "dust_w must be a positive integer.");
            }
            break;
        case 0x800:
            arguments->dust_t = atoi(arg);
            if (arguments->dust_t <= 0) {
                argp_error(state, "dust_t must be a positive integer.");
            }
            break;
        case 'q':
            arguments->min_qscore = (float)atof(arg);
            break;
        case 'x':
            arguments->recurse = -1;  // 0: stops recursion
            break;
        case 'H':
            arguments->reheader = 1;
            break;
        case 'B':
            arguments->write_bam = 1;
            break;
        case 'v':
            arguments->verbose = 1;
            break;
        case 't':
            arguments->threads = atoi(arg);
            break;
        case 'e':
            arguments->force_error = 1;
            break; 
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            arguments->files = &state->argv[state->next - 1];
            state->next = state->argc;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};


arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    args.perread = NULL;
    args.perfile = NULL;
    args.runids = NULL;
    args.basecallers = NULL;
    args.sample = "";
    args.min_length = 0;
    args.max_length = (size_t)-1;
    args.min_qscore = 0;
    args.recurse = 1; // always allow descent into TLD
    args.demultiplex_dir = NULL;
    args.histograms = "fastcat-histograms";
    args.dust = 0;
    args.max_dust = 0.95;
    args.dust_w = 64;
    args.dust_t = 20;
    args.reheader = 0;
    args.write_bam = 0;
    args.threads = 1;
    args.reads_per_file = 0;
    args.force_error = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}
