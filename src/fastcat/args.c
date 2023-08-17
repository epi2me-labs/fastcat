#include <stdlib.h>
#include <argp.h>

#include "args.h"
#include "../version.h"

const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
"fastcat -- concatenate and summarise .fastq(.gz) files.\
\vInput files may be given on stdin by specifing the input as '-'. \
Also accepts directories as input and looks for .fastq(.gz) files in \
the top-level directory. Recurses into sub-directories when the \
-x option is given. The command \
will exit non-zero if any file encountered cannot be read.";
static char args_doc[] = "reads1.fastq(.gz) reads2.fastq(.gz) dir-with-fastq ...";
static struct argp_option options[] = {
    {"read", 'r', "READ SUMMARY",  0,
        "Per-read summary output", 0},
    {"file", 'f', "FILE SUMMARY",  0,
        "Per-file summary output", 0},
    {"sample", 's', "SAMPLE NAME",   0,
        "Sample name (if given, adds a 'sample_name' column).", 0},
    {"demultiplex", 'd', "OUT DIR",  0,
        "Separate barcoded samples using fastq header information. Option value is top-level output directory.", 0},
    {"min_length", 'a', "MIN READ LENGTH", 0,
        "minimum read length to output (excluded reads remain listed in summaries).", 0},
    {"max_length", 'b', "MAX READ LENGTH", 0,
        "maximum read length to output (excluded reads remain listed in summaries).", 0},
    {"min_qscore", 'q', "MIN READ QSCOROE", 0,
        "minimum read Qscore to output (excluded reads remain listed in summaries).", 0},
    {"recurse", 'x', 0, 0,
        "Search directories recursively for '.fastq', '.fq', '.fastq.gz', and '.fq.gz' files.", 0},
    {"reheader", 'H', 0, 0,
        "Rewrite fastq header comments as SAM tags (useful for passing through minimap2).", 0},
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
        case 's':
            arguments->sample = arg;
            break;
        case 'a':
            arguments->min_length = atoi(arg);
            break;
        case 'b':
            arguments->max_length = atoi(arg);
            break;
        case 'd':
            arguments->demultiplex_dir = arg;
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
    args.sample = "";
    args.min_length = 0;
    args.max_length = (size_t)-1;
    args.min_qscore = 0;
    args.recurse = 1; // always allow descent into TLD
    args.demultiplex_dir = NULL;
    args.reheader = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}
