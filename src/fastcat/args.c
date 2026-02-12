#include <stdlib.h>
#include <string.h>
#include <argp.h>

#include "args.h"
#include "../version.h"
#include "../fc_log.h"
#include "../fc_input_list.h"

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
    {"output", 'o', "DIRECTORY", 0,
        "Output root directory (default: fastcat-qc). Must not already exist.", 0},
    {"sample", 's', "SAMPLE NAME",   0,
        "Sample name (if given, adds a 'sample_name' column).", 0},
    {"reads-per-file", 'c', "NUM", 0,
        "Split reads into files with a set number of reads (default: single file).", 0},
    {"reheader", 'H', 0, 0,
        "Rewrite fastq header comments as SAM tags (useful for passing through minimap2).", 0},
    {"bam_out", 'B', 0, 0,
        "Output data as unaligned BAM.", 0},
    {"verbose", 'v', 0, 0,
        "Verbose output.", 0},

    {0, 0, 0, 0,
        "Output file selection:", 0},
    {"per-read", 'p', 0, 0,
        "Write per-read TSV (large). Disabled by default.", 0},
    {"demultiplex", 'd', 0, 0,
        "Separate barcoded samples using fastq header information; outputs go under <output>/barcodeXXXX/.", 0},

    {0, 0, 0, 0,
        "Read filtering options:", 0},
    {"min-length", 'a', "MIN READ LENGTH", 0,
        "minimum read length to output (excluded reads remain listed in summaries).", 0},
    {"max-length", 'b', "MAX READ LENGTH", 0,
        "maximum read length to output (excluded reads remain listed in summaries).", 0},
    {"min-qscore", 'q', "MIN READ QSCORE", 0,
        "minimum read Qscore to output (excluded reads remain listed in summaries).", 0},
    {"dust", 0x500, 0, 0,
        "Enable DUST filtering of reads (default: disabled).", 0},

    { 0, 0, 0, 0,
        "Advanced cleaning options:", 0},
    {"max-dust", 0x600, "MAX DUST", 0,
        "Maximum proportion of low-complexity regions to allow in reads (default: 0.95).", 0},
    {"dust-w", 0x700, "DUST W", 0,
        "Window size for DUST filtering (default: 64).", 0},
    {"dust-t", 0x800, "DUST T", 0,
        "Threshold for DUST filtering (default: 20).", 0},

    {0, 0, 0, 0,
        "Logging options:", 0},
    {"log-level", 0x900, "LEVEL", 0,
        "Log level: error, warn, info, debug (default: info).", 0},
    {"log-file",  0x901, "PATH", 0,
        "Write logs to PATH instead of stderr.", 0},
    {"log-mono",  0x902, 0, 0,
        "Force monochrome logging output (disable ANSI colors).", 0},
    { 0 }
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    switch (key) {
        case 's':
            arguments->sample = arg;
            break;
        case 'o':
            arguments->output_dir = arg;
            break;
        case 'p':
            arguments->per_read = true;
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
            arguments->demultiplex = true;
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
        case 0x900:
            arguments->log_level = fc_log_level_from_string(arg, FC_LOG_INFO);
            if (arguments->log_level == FC_LOG_INFO && arg && strcasecmp(arg, "info") != 0) {
                argp_error(state, "Unknown log level '%s' (expected error|warn|info|debug).", arg);
            }
            break;
        case 0x901:
            arguments->log_file = arg;
            break;
        case 0x902:
            arguments->log_mono = true;
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
    args.output_dir = "fastcat-qc";
    args.sample = "";
    args.min_length = 0;
    args.max_length = (size_t)-1;
    args.min_qscore = 0;
    args.recurse = 1; // always allow descent into TLD
    args.demultiplex = false;
    args.per_read = false;
    args.dust = 0;
    args.max_dust = 0.95;
    args.dust_w = 64;
    args.dust_t = 20;
    args.reheader = 0;
    args.write_bam = 0;
    args.threads = 1;
    args.reads_per_file = 0;
    args.force_error = 0;
    args.log_level = FC_LOG_INFO;
    args.log_file = NULL;
    args.log_mono = false;
    args.verbose = false;
    args.files_from_stdin = false;
    args.n_files = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);

    // If files are provided via stdin ("-"), read them now.
    if (args.files && args.files[0] && strcmp(args.files[0], "-") == 0) {
        size_t n = 0;
        char **paths = fc_read_paths_from_stdin(&n);
        args.files_from_stdin = true;
        args.n_files = n;
        args.files = paths;
    } else if (args.files) {
        while (args.files[args.n_files]) args.n_files++;
    }
    return args;
}

void destroy_args(arguments_t *args) {
    if (!args) return;
    if (args->files_from_stdin && args->files) {
        fc_free_path_list(args->files, args->n_files);
        args->files = NULL;
    }
}
