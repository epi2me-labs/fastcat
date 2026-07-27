#include "cli/catstats/args.h"

#include <argp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "core/fc_log.h"
#include "core/fc_time_bin.h"
#include "version.h"

const char *argp_program_bug_address = "support@nanoporetech.com";

enum {
    KEY_DUST = 0x1000,
    KEY_MAX_DUST,
    KEY_DUST_W,
    KEY_DUST_T,
    KEY_FOFN,
    KEY_INPUT_FORMAT,
    KEY_RECALC_QUAL,
    KEY_TAG_NAME,
    KEY_TAG_VALUE,
    KEY_HAPLOTYPE,
    KEY_COVERAGE_BEDS,
    KEY_COVERAGE_NAMES,
    KEY_THRESHOLDS,
    KEY_SEGMENTS,
    KEY_COVERAGE,
    KEY_POLYA,
    KEY_POLYA_COVER,
    KEY_POLYA_QUAL,
    KEY_POLYA_REV,
    KEY_LOG_LEVEL,
    KEY_LOG_FILE,
    KEY_LOG_MONO,
    KEY_MAX_LENGTH_LONG,
    KEY_LENGTH_QSCORE_K,
    KEY_LENGTH_QSCORE_QGRID,
    KEY_LENGTH_QSCORE_QSMOOTH,
    KEY_TIME_BIN_SIZE,
    KEY_TIME_BIN_K,
    KEY_HTML_REPORT,
    KEY_NO_LEGACY_OUTPUTS
};

typedef struct {
    fc_program_kind program;
    fc_catstats_args *args;
    int tag_items;
    bool tag_given;
    bool hp_given;
} parse_context;

// Shared option macros to avoid duplication between fastcat and bamstats
#define SHARED_GENERAL_OPTIONS \
    { "force-error", 'e', 0, 0, \
        "Exit with non-zero status if any files, or records, contained errors.", 0 }, \
        { "verbose", 'v', 0, 0, "Verbose output.", 0 }

#define SHARED_OUTPUT_OPTIONS \
    { "output", 'o', "DIRECTORY", 0, \
        "Output root directory (default: fastcat-qc). Must not already exist.", 10 }, \
        { "demultiplex", 'd', 0, 0, \
            "Separate barcoded samples using fastq header information; outputs go under <output>/barcodeXXXX/.", 10 }, \
        { "sample", 's', "SAMPLE NAME", 0, \
            "Sample identifier (defaults to 'SAMPLE' if files don't provide sample_id metadata).", 10 }, \
        { "per-read", 'p', 0, 0, "Write per-read TSV (large). Disabled by default.", 12 }, \
        { "html-report", KEY_HTML_REPORT, 0, 0, \
            "Write a single-sample HTML report using aggregate fastcat histograms.", 13 }, \
        { "no-legacy-outputs", KEY_NO_LEGACY_OUTPUTS, 0, 0, \
            "Disable legacy counter and histogram summary outputs; summary.json is still written.", 13 }

#define SHARED_FILTER_OPTIONS \
    { "min-length", 'a', "MIN READ LENGTH", 0, \
        "Minimum read length to output (excluded reads remain listed in summaries).", 20 }, \
        { "min-qscore", 'q', "MIN READ QSCORE", 0, \
            "Minimum read Qscore to output (excluded reads remain listed in summaries).", 20 }, \
        { "dust", KEY_DUST, 0, 0, "Enable DUST filtering of reads (default: disabled).", 20 }, \
        { "max-dust", KEY_MAX_DUST, "MAX DUST", 0, \
            "Maximum proportion of low-complexity regions to allow in reads (default: 0.95).", 20 }, \
        { "dust-w", KEY_DUST_W, "DUST W", 0, "Window size for DUST filtering (default: 64).", 20 }, \
        { "dust-t", KEY_DUST_T, "DUST T", 0, "Threshold for DUST filtering (default: 20).", 20 }, \
        { "read-group", 'g', "RG", 0, "Only process reads from given read group.", 20 }, \
        { "tag-name", KEY_TAG_NAME, "TN", 0, "Only process reads with a given tag (see --tag-value).", 20 }, \
        { "tag-value", KEY_TAG_VALUE, "VAL", 0, "Only process reads with a given tag value.", 20 }, \
        { "haplotype", KEY_HAPLOTYPE, "VAL", 0, \
            "Only process reads from a given haplotype. Equivalent to --tag-name HP --tag-value VAL.", 20 }, \
        { "recalc-qual", KEY_RECALC_QUAL, 0, 0, \
            "Force recomputing mean quality, else use 'qs' tag if present.", 20 }

#define SHARED_HISTOGRAM_OPTIONS \
    { "length-qscore-k", KEY_LENGTH_QSCORE_K, "K", 0, \
        "Sketch size for length modeling per qscore bin (higher => more memory, less sketch noise; default: 400).", 30 }, \
        { "length-qscore-qgrid", KEY_LENGTH_QSCORE_QGRID, "N", 0, \
            "Quantile samples converted into each density curve (higher => finer detail, more CPU; default: 257).", 30 }, \
        { "length-qscore-qsmooth", KEY_LENGTH_QSCORE_QSMOOTH, "W", 0, \
            "Odd moving-average window over quantile slopes (higher => smoother/less spiky, but flatter peaks; default: 5).", 30 }

#define SHARED_TIME_OPTIONS \
    { "time-bin-size", KEY_TIME_BIN_SIZE, "SECONDS", 0, \
        "Width of each time bin in seconds (higher => coarser time plots, lower => more rows; default: 3600).", 35 }, \
        { "time-bin-k", KEY_TIME_BIN_K, "K", 0, \
            "Sketch size for per-bin length/qscore quantiles (higher => more memory, less sketch noise; default: 200).", 35 }

#define SHARED_LOGGING_OPTIONS \
    { "log-level", KEY_LOG_LEVEL, "LEVEL", 0, \
        "Log level: error, warn, info, debug (default: info).", 40 }, \
        { "log-file", KEY_LOG_FILE, "PATH", 0, "Write logs to PATH instead of stderr.", 40 }, \
        { "log-mono", KEY_LOG_MONO, 0, 0, "Force monochrome logging output (disable ANSI colors).", 40 }

static char fastcat_doc[] = "fastcat -- concatenate and summarise FASTQ or unaligned BAM reads.\
\vReads input data from stdin when '-' is provided as an input path. \
Use --input-fmt to choose FASTQ (default) or unaligned BAM input. \
Use --fofn to interpret each input path as a file-of-filenames source. \
Also accepts directories as input and looks for files matching the selected \
input format in the top-level directory. Recurses into sub-directories when the \
-x option is given. The command \
will exit non-zero if any file encountered cannot be read.";

static char fastcat_args_doc[] = "reads1 reads2 dir-with-reads ...";

static struct argp_option fastcat_options[] = {
    { 0, 0, 0, 0, "General options:", 0 },
    { "recurse", 'x', 0, 0,
        "Search directories recursively for files matching the selected input format.", 0 },
    { "fofn", KEY_FOFN, 0, 0,
        "Treat each input as a file-of-filenames source (use '-' to read that list from stdin).", 0 },
    { "input-fmt", KEY_INPUT_FORMAT, "fastq|bam", 0,
        "Input format for fastq mode (default: fastq). Use 'bam' for unaligned BAM input.", 0 },
    { "threads", 't', "THREADS", 0,
        "Number of threads for output compression (only with --bam-out).", 0 },
    SHARED_GENERAL_OPTIONS,

    { 0, 0, 0, 0, "Output options:", 10 },
    SHARED_OUTPUT_OPTIONS,
    { "bam-out", 'B', 0, 0, "Output data as unaligned BAM.", 11 },
    { "reads-per-file", 'c', "NUM", 0,
        "Split reads into files with a set number of reads (default: single file).", 11 },

    { 0, 0, 0, 0, "Read filtering options:", 20 },
    SHARED_FILTER_OPTIONS,
    { "max-length", 'b', "MAX READ LENGTH", 0,
        "Maximum read length to output (excluded reads remain listed in summaries).", 20 },

    { 0, 0, 0, 0, "Histogram Sketch Options: controls summary.json heatmap detail vs smoothness.", 30 },
    SHARED_HISTOGRAM_OPTIONS,

    { 0, 0, 0, 0, "Time Summary Options: controls summary.json time-series resolution and quantile accuracy.", 35 },
    SHARED_TIME_OPTIONS,

    { 0, 0, 0, 0, "Logging options:", 40 },
    SHARED_LOGGING_OPTIONS,
    { 0 }
};

static char bam_doc[] = "fastcat bam -- summarise reads/alignments in an input BAM file.\
\vCreates summary files containing statistics for reads and alignments in \
the input BAM file.";

static char bam_args_doc[] = "<reads.bam>";

static struct argp_option bam_options[] = {
    { 0, 0, 0, 0, "General options:", 0 },
    { "region", 'r', "chr:start-end", 0, "Genomic region to process.", 0 },
    { "bed", 'b', "BEDFILE", 0, "BED file for regions to process.", 0 },
    { "threads", 't', "THREADS", 0,
        "Number of threads (for BAM decompression and BED output compression).", 0 },
    SHARED_GENERAL_OPTIONS,

    { 0, 0, 0, 0, "Output options:", 10 },
    SHARED_OUTPUT_OPTIONS,

    { 0, 0, 0, 0, "Read filtering options:", 20 },
    SHARED_FILTER_OPTIONS,
    { "max-length", KEY_MAX_LENGTH_LONG, "MAX READ LENGTH", 0,
        "Maximum read length to output (excluded reads remain listed in summaries).", 20 },

    { 0, 0, 0, 0, "Histogram Sketch Options: controls summary.json heatmap detail vs smoothness.", 30 },
    SHARED_HISTOGRAM_OPTIONS,

    { 0, 0, 0, 0, "Time Summary Options: controls summary.json time-series resolution and quantile accuracy.", 35 },
    SHARED_TIME_OPTIONS,

    { 0, 0, 0, 0, "Coverage calculation options:", 40 },
    { "coverage", KEY_COVERAGE, 0, 0, "Enable coverage calculations.", 40 },
    { "coverage-beds", KEY_COVERAGE_BEDS, "BEDFILE ...", 0,
        "BED file(s) for calculating coverage (space-separated list).", 40 },
    { "coverage-names", KEY_COVERAGE_NAMES, "NAME ...", 0,
        "Name(s) for coverage BED file(s) (space-separated list).", 40 },
    { "thresholds", KEY_THRESHOLDS, "VALUES ...", 0,
        "Coverage thresholds to produce sparse cumulative distribution (space-separated integers).", 40 },
    { "segments", KEY_SEGMENTS, "LENGTH ...", 0,
        "Segment length(s) for which to produce outputs (space-separated integers).", 40 },

    { 0, 0, 0, 0, "Poly-A options:", 50 },
    { "poly-a", KEY_POLYA, 0, 0, "Enable poly-A tail length histogram.", 50 },
    { "poly-a-cover", KEY_POLYA_COVER, "PCT_COVERAGE", 0,
        "Reference alignment coverage for acceptance of read (default: 95).", 50 },
    { "poly-a-qual", KEY_POLYA_QUAL, "QUAL", 0,
        "Read mean Q score for acceptance of read (default: 10).", 50 },
    { "poly-a-rev", KEY_POLYA_REV, 0, 0,
        "Allow reverse alignments (useful for cDNA; default is direct RNA).", 50 },

    { 0, 0, 0, 0, "Logging options:", 60 },
    SHARED_LOGGING_OPTIONS,
    { 0 }
};

static fc_catstats_args default_args(void)
{
    fc_catstats_args args = { 0 };
    args.output_dir = "fastcat-qc";
    args.sample = NULL;
    args.demultiplex = false;
    args.per_read = false;
    args.html_report = false;
    args.no_legacy_outputs = false;
    args.threads = 1;
    args.verbose = false;
    args.force_error = false;
    args.log_level = FC_LOG_INFO;
    args.log_file = NULL;
    args.log_mono = false;
    args.length_qscore_k = 400;
    args.length_qscore_qgrid = 257;
    args.length_qscore_qsmooth = 5;
    args.time_bin_size_s = FC_TIME_BIN_DEFAULT_SIZE_S;
    args.time_bin_k = FC_TIME_BIN_DEFAULT_K;

    args.min_length = 0;
    args.max_length = (size_t)-1;
    args.min_qscore = 0.0f;
    args.dust = false;
    args.max_dust = 0.95;
    args.dust_w = 64;
    args.dust_t = 20;
    args.read_group = NULL;
    args.tag_name[0] = '\0';
    args.tag_name[1] = '\0';
    args.tag_value = -1;
    args.force_recalc_qual = false;

    args.recurse = 1;
    args.write_bam = 0;
    args.fofn = false;
    args.input_kind = FC_FILE_FASTQ;
    args.reads_per_file = 0;
    args.input_paths.paths = NULL;
    args.input_paths.n_paths = 0;

    args.bam = NULL;
    args.region = NULL;
    args.bed = NULL;

    args.coverage = false;
    args.coverage_beds = NULL;
    args.n_coverage_beds = 0;
    args.coverage_names = NULL;
    args.n_coverage_names = 0;
    args.coverage_thresholds = NULL;
    args.n_coverage_thresholds = 0;
    args.segments = NULL;
    args.n_segments = 0;

    args.poly_a = false;
    args.poly_a_cover = 95.0f;
    args.poly_a_qual = 10.0f;
    args.poly_a_rev = false;
    return args;
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    parse_context *ctx = (parse_context *)state->input;
    fc_catstats_args *arguments = ctx->args;
    size_t size_value = 0;
    int int_value = 0;
    float float_value = 0.0f;
    double double_value = 0.0;

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
    case KEY_HTML_REPORT:
        arguments->html_report = true;
        break;
    case KEY_NO_LEGACY_OUTPUTS:
        arguments->no_legacy_outputs = true;
        break;
    case 'a':
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --min-length value: '%s'.", arg);
        }
        arguments->min_length = size_value;
        break;
    case 'b':
        if (ctx->program == FC_PROGRAM_FASTCAT) {
            if (!fc_parse_size_strict(arg, &size_value)) {
                argp_error(state, "Invalid --max-length value: '%s'.", arg);
            }
            arguments->max_length = size_value;
        } else {
            arguments->bed = arg;
        }
        break;
    case KEY_MAX_LENGTH_LONG:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --max-length value: '%s'.", arg);
        }
        arguments->max_length = size_value;
        break;
    case KEY_LENGTH_QSCORE_K:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --length-qscore-k value: '%s'.", arg);
        }
        if (size_value < 8) {
            argp_error(state, "--length-qscore-k must be >= 8.");
        }
        arguments->length_qscore_k = size_value;
        break;
    case KEY_LENGTH_QSCORE_QGRID:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --length-qscore-qgrid value: '%s'.", arg);
        }
        if (size_value < 2) {
            argp_error(state, "--length-qscore-qgrid must be >= 2.");
        }
        arguments->length_qscore_qgrid = size_value;
        break;
    case KEY_LENGTH_QSCORE_QSMOOTH:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --length-qscore-qsmooth value: '%s'.", arg);
        }
        if (size_value == 0 || (size_value % 2) == 0) {
            argp_error(state, "--length-qscore-qsmooth must be odd and >= 1.");
        }
        arguments->length_qscore_qsmooth = size_value;
        break;
    case KEY_TIME_BIN_SIZE:
        if (!fc_parse_double_strict(arg, &double_value)) {
            argp_error(state, "Invalid --time-bin-size value: '%s'.", arg);
        }
        if (!isfinite(double_value) || double_value <= 0.0) {
            argp_error(state, "--time-bin-size must be positive.");
        }
        arguments->time_bin_size_s = double_value;
        break;
    case KEY_TIME_BIN_K:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --time-bin-k value: '%s'.", arg);
        }
        if (size_value < 8) {
            argp_error(state, "--time-bin-k must be >= 8.");
        }
        arguments->time_bin_k = size_value;
        break;
    case 'c':
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --reads-per-file value: '%s'.", arg);
        }
        arguments->reads_per_file = size_value;
        break;
    case 'd':
        arguments->demultiplex = true;
        break;
    case KEY_DUST:
        arguments->dust = true;
        break;
    case KEY_MAX_DUST:
        if (!fc_parse_double_strict(arg, &double_value)) {
            argp_error(state, "Invalid --max-dust value: '%s'.", arg);
        }
        arguments->max_dust = double_value;
        if (arguments->max_dust < 0.0 || arguments->max_dust > 1.0) {
            argp_error(state, "max_dust must be between 0 and 1.");
        }
        break;
    case KEY_DUST_W:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --dust-w value: '%s'.", arg);
        }
        arguments->dust_w = size_value;
        if (arguments->dust_w == 0) {
            argp_error(state, "dust_w must be a positive integer.");
        }
        break;
    case KEY_DUST_T:
        if (!fc_parse_size_strict(arg, &size_value)) {
            argp_error(state, "Invalid --dust-t value: '%s'.", arg);
        }
        arguments->dust_t = size_value;
        if (arguments->dust_t == 0) {
            argp_error(state, "dust_t must be a positive integer.");
        }
        break;
    case KEY_LOG_LEVEL:
        arguments->log_level = fc_log_level_from_string(arg, FC_LOG_INFO);
        if (arguments->log_level == FC_LOG_INFO && arg && strcasecmp(arg, "info") != 0) {
            argp_error(state, "Unknown log level '%s' (expected error|warn|info|debug).", arg);
        }
        break;
    case KEY_LOG_FILE:
        arguments->log_file = arg;
        break;
    case KEY_LOG_MONO:
        arguments->log_mono = true;
        break;
    case 'q':
        if (!fc_parse_float_strict(arg, &float_value)) {
            argp_error(state, "Invalid --min-qscore value: '%s'.", arg);
        }
        arguments->min_qscore = float_value;
        break;
    case 'x':
        arguments->recurse = -1; // 0 stops recursion
        break;
    case 'B':
        arguments->write_bam = 1;
        break;
    case 'v':
        arguments->verbose = true;
        break;
    case 't':
        if (!fc_parse_int_strict(arg, &int_value)) {
            argp_error(state, "Invalid --threads value: '%s'.", arg);
        }
        if (int_value <= 0) {
            argp_error(state, "--threads must be > 0.");
        }
        arguments->threads = int_value;
        break;
    case 'e':
        arguments->force_error = true;
        break;
    case KEY_FOFN:
        arguments->fofn = true;
        break;
    case KEY_INPUT_FORMAT:
        if (strcasecmp(arg, "fastq") == 0) {
            arguments->input_kind = FC_FILE_FASTQ;
        } else if (strcasecmp(arg, "bam") == 0) {
            arguments->input_kind = FC_FILE_BAM;
        } else {
            argp_error(state, "Unsupported --input-fmt '%s' (expected fastq or bam).", arg);
        }
        break;
    case KEY_RECALC_QUAL:
        arguments->force_recalc_qual = true;
        break;
    case 'g':
        arguments->read_group = arg;
        break;
    case KEY_TAG_NAME:
        if (strlen(arg) != 2) {
            argp_error(state, "Tag name should be a two-letter code, received: '%s'.", arg);
        }
        memcpy(arguments->tag_name, arg, 2 * sizeof(char));
        ctx->tag_items += 1;
        ctx->tag_given = true;
        break;
    case KEY_TAG_VALUE:
        if (!fc_parse_int_strict(arg, &int_value)) {
            argp_error(state, "Invalid --tag-value value: '%s'.", arg);
        }
        arguments->tag_value = int_value;
        ctx->tag_items += 1;
        ctx->tag_given = true;
        break;
    case KEY_HAPLOTYPE:
        memcpy(arguments->tag_name, "HP", 2 * sizeof(char));
        if (!fc_parse_int_strict(arg, &int_value)) {
            argp_error(state, "Invalid --haplotype value: '%s'.", arg);
        }
        arguments->tag_value = int_value;
        ctx->tag_items += 2;
        ctx->hp_given = true;
        break;
    case 'r':
        arguments->region = arg;
        break;
    case KEY_COVERAGE_BEDS:
        slurp_args((const char ***)&arguments->coverage_beds, &arguments->n_coverage_beds, arg, state);
        break;
    case KEY_COVERAGE_NAMES:
        slurp_args((const char ***)&arguments->coverage_names, &arguments->n_coverage_names, arg, state);
        break;
    case KEY_THRESHOLDS:
        slurp_ints(&arguments->coverage_thresholds, &arguments->n_coverage_thresholds, arg, state);
        break;
    case KEY_SEGMENTS:
        slurp_ints(&arguments->segments, &arguments->n_segments, arg, state);
        break;
    case KEY_COVERAGE:
        arguments->coverage = true;
        break;
    case KEY_POLYA:
        arguments->poly_a = true;
        break;
    case KEY_POLYA_COVER:
        if (!fc_parse_float_strict(arg, &float_value)) {
            argp_error(state, "Invalid --poly-a-cover value: '%s'.", arg);
        }
        if (float_value < 0.0 || float_value > 100.0) {
            argp_error(state, "--poly-a-cover must be between 0.0 and 100.0 (percentage).");
        }
        arguments->poly_a_cover = float_value;
        break;
    case KEY_POLYA_QUAL:
        if (!fc_parse_float_strict(arg, &float_value)) {
            argp_error(state, "Invalid --poly-a-qual value: '%s'.", arg);
        }
        if (float_value < 0.0 || float_value > 60.0) {
            argp_error(state, "--poly-a-qual must be between 0.0 and 60.0.");
        }
        arguments->poly_a_qual = float_value;
        break;
    case KEY_POLYA_REV:
        arguments->poly_a_rev = true;
        break;
    case ARGP_KEY_ARG:
        if (ctx->program == FC_PROGRAM_FASTCAT) {
            arguments->input_paths.paths = &state->argv[state->next - 1];
            state->next = state->argc;
        } else {
            if (state->arg_num == 0) {
                arguments->bam = arg;
            } else {
                argp_usage(state);
            }
        }
        break;
    case ARGP_KEY_NO_ARGS:
        if (ctx->program == FC_PROGRAM_BAMSTATS || ctx->program == FC_PROGRAM_FASTCAT) {
            argp_usage(state);
        }
        break;
    case ARGP_KEY_END:
        if (ctx->tag_items % 2 > 0) {
            argp_error(state, "Both or neither of --tag-name and --tag-value must be given.");
        }
        if (ctx->tag_given && ctx->hp_given) {
            argp_error(
                state, "If --haplotype is given neither --tag-name nor --tag-value should be provided.");
        }
        if (ctx->program == FC_PROGRAM_FASTCAT) {
            if (!arguments->input_paths.paths) {
                argp_error(state, "Missing input path(s). Use '-' to read data from stdin.");
                return 0; // argp_error doesn't return, but helps static analyzer
            }
            while (arguments->input_paths.paths[arguments->input_paths.n_paths]) {
                arguments->input_paths.n_paths++;
            }
        } else {
            if (!arguments->bam) {
                argp_error(state, "Missing <reads.bam>");
            }
            if (arguments->region && arguments->bed) {
                argp_error(state, "--region and --bed are mutually exclusive.");
            }
            if (arguments->n_coverage_beds && arguments->n_coverage_beds != arguments->n_coverage_names) {
                argp_error(
                    state,
                    "Mismatched counts: --coverage-beds (%zu) vs --coverage-names (%zu).",
                    arguments->n_coverage_beds,
                    arguments->n_coverage_names);
            }
        }
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp fastcat_argp = { fastcat_options, parse_opt, fastcat_args_doc, fastcat_doc, 0, 0, 0 };
static struct argp bam_argp = { bam_options, parse_opt, bam_args_doc, bam_doc, 0, 0, 0 };

fc_catstats_args fc_parse_catstats_args(fc_program_kind program, int argc, char **argv)
{
    fc_catstats_args args = default_args();
    parse_context ctx = {
        .program = program,
        .args = &args,
        .tag_items = 0,
        .tag_given = false,
        .hp_given = false
    };
    struct argp *tool_argp = (program == FC_PROGRAM_FASTCAT) ? &fastcat_argp : &bam_argp;
    argp_parse(tool_argp, argc, argv, 0, 0, &ctx);
    return args;
}

void fc_destroy_catstats_args(fc_catstats_args *args)
{
    if (!args) {
        return;
    }
    free(args->coverage_beds);
    free(args->coverage_names);
    free(args->coverage_thresholds);
    free(args->segments);
}
