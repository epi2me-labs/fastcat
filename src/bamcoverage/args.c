#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>

#include "args.h"
#include "common.h"
#include "../version.h"

const char *argp_program_bug_address = "support@nanoporetech.com";

static char doc[] =
"bamcoverage -- summarise coverage statistics from a BAM file.\n"
"\vThe program creates coverage traces in BED files, and various summary files.";

static char args_doc[] = "<reads.bam>";

static struct argp_option options[] = {
    { 0, 0, 0, 0, "General options:", 0 },
    { "beds", 'b', "BEDFILE ...", 0, "BED files for regions (space-separated list).", 0 },
    { "names", 'n', "NAME...", 0, "Names associated with the BED files (space-separated list).", 0 },
    { "threads", 't', "THREADS", 0, "Number of threads for BAM processing.", 0 },
    { 0 }
};


static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    arguments_t *a = state->input;

    switch (key) {
        case 'b':
            slurp_args(&a->beds, &a->n_beds, arg, state);
            break;
        case 'n':
            slurp_args(&a->bed_names, &a->n_bed_names, arg, state);
            break;
        case 't':
            a->threads = atoi(arg);
            if (a->threads <= 0) argp_error(state, "THREADS must be > 0");
            break;

        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                a->bam = arg;
                break;
            }
            argp_usage(state);
            break;

        case ARGP_KEY_NO_ARGS:
            argp_usage(state);
            break;

        case ARGP_KEY_END:
            if (!a->bam) argp_error(state, "Missing <reads.bam>");
            if (a->n_beds && a->n_bed_names != a->n_beds) {
                argp_error(state, "Mismatched counts: --bed (%zu) vs --bed-name (%zu).",
                           a->n_beds, a->n_bed_names);
            }
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


arguments_t parse_arguments(int argc, char **argv) {
    arguments_t a = {
        .bam = NULL,
        .beds = NULL, .n_beds = 0,
        .bed_names = NULL, .n_bed_names = 0,
        .threads = 1,
    };
    argp_parse(&argp, argc, argv, 0, 0, &a);
    return a;
}
