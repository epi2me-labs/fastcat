#include <stdlib.h>
#include <string.h>
#include <argp.h>

#include "args.h"
#include "../version.h"

const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
"fastlint -- apply sdust algorithm to input files.\
\vThe program removes low-complexity reads from the input stream.\n";
static char args_doc[] = "<reads.fastq>";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:", 0},
    {"threshold", 't', "THRESHOLD", 0,
        "Threshold for repetition (default: 20).", 0},
    {"window", 'w', "WINDOW", 0,
        "Window size (default: 64).", 0},
    {"max_proportion", 'p', "PROPORTION", 0,
        "Maximum allowable proportion of masked bases in a read to keep the read (default: 0.95).", 0},
    { 0 }
};


static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 't':
            arguments->t = atoi(arg);
            break;
        case 'w':
            arguments->w = atoi(arg);
            break;
        case 'p':
            arguments->max_proportion = atof(arg);
            if (arguments->max_proportion < 0.0 || arguments->max_proportion > 1.0) {
                argp_error(state, "Proportion must be between 0.0 and 1.0.");
            }
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            // Count positional arguments
            if (arguments->nfiles >= MAX_INPUT_FILES) {
                argp_error(state, "Too many input files (max %d)", MAX_INPUT_FILES);
            }
            arguments->fastq[arguments->nfiles++] = arg;
            break;
        case ARGP_KEY_END:
            if (arguments->nfiles == 0)
            argp_usage(state);
            arguments->fastq[arguments->nfiles] = NULL; // NULL-terminate
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    memset(&args, 0, sizeof(arguments_t)); // initialize fastqs
    args.t = 20;
    args.w = 64;
    args.max_proportion = 0.95;
    args.nfiles = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}
