// bamindex build program

#include <argp.h>
#include <err.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "../common.h"
#include "../fc_log.h"
#include "../version.h"
#include "index.h"

typedef struct arguments {
    const char *bam;
    int threads;
    int chunk_size;
} arguments_t;

static char doc[] = "bamindex build -- create a BAM index corresponding to batches of records.\
\vThe program creates a simple index of file offsets for each of every \
(n * M)th alignment record. No care is taken to keep records corresponding \
to the same query together, or any other such niceities. Its intended to \
be used simply with unaligned, unsorted BAMs.";
static char args_doc[] = "<reads.bam>";
static struct argp_option options[] = {
    { 0, 0, 0, 0,
        "General options:", 0 },
    { "threads", 't', "THREADS", 0,
        "Number of threads for BAM processing.", 0 },
    { "chunk-size", 'c', "SIZE", 0,
        "Number of records in a chunk.", 0 },
    { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    arguments_t *arguments = state->input;
    int int_value = 0;
    switch (key) {
    case 't':
        if (!fc_parse_int_strict(arg, &int_value)) {
            argp_error(state, "Invalid --threads value: '%s'.", arg);
        }
        if (int_value <= 0) {
            argp_error(state, "--threads must be > 0.");
        }
        arguments->threads = int_value;
        break;
    case 'c':
        if (!fc_parse_int_strict(arg, &int_value)) {
            argp_error(state, "Invalid --chunk-size value: '%s'.", arg);
        }
        if (int_value <= 0) {
            argp_error(state, "--chunk-size must be > 0.");
        }
        arguments->chunk_size = int_value;
        break;
    case ARGP_KEY_NO_ARGS:
        argp_usage(state);
        break;
    case ARGP_KEY_ARG:
        if (state->arg_num == 0) {
            arguments->bam = arg;
            break;
        }
        break;
    case ARGP_KEY_END:
        if (state->arg_num != 1) {
            argp_usage(state);
        }
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

static arguments_t parse_arguments(int argc, char **argv)
{
    arguments_t args;
    args.bam = NULL;
    args.threads = 1;
    args.chunk_size = 1;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}

void index_build(const char *filename, const char *output_fname, int threads, size_t every)
{
    htsFile *fp = hts_open(filename, "r");
    bam_hdr_t *h = sam_hdr_read(fp);
    if (fp == NULL || h == NULL) {
        FC_LOG_FATAL("Could not open %s", filename);
    }

    FILE *out_fp = fopen(output_fname, "wb");
    if (!out_fp) {
        FC_LOG_ERROR("Cannot open output file: %s.", output_fname);
        bam_hdr_destroy(h);
        hts_close(fp);
        exit(EXIT_FAILURE);
    }
    bc_idx_t *idx = bc_idx_init1(every);
    size_t rtn;
    if ((rtn = bc_idx_write_header(out_fp, idx)) > 0) {
        FC_LOG_ERROR("Failed to write header to index. Error %zu.", rtn);
    }

    htsThreadPool p = { NULL, 0 };
    if (threads > 1) {
        p.pool = hts_tpool_init(threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    int ret = 0;
    int i = 0;
    bam1_t *b = bam_init1();
    size_t file_offset = bgzf_tell(fp->fp.bgzf);
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        if ((i % every) != 0) {
            file_offset = bgzf_tell(fp->fp.bgzf);
            i++;
            continue;
        }
        if (i % 100000 == 0) {
            FC_LOG_INFO("Record %d %zu", i, file_offset);
        }
        if (bc_idx_write(out_fp, idx, file_offset, bam_get_qname(b)) < 0) {
            FC_LOG_FATAL("Failed to write records to index.");
        }
        file_offset = bgzf_tell(fp->fp.bgzf);
        i++;
    }

    if (ret < -1) {
        FC_LOG_FATAL("Error reading BAM file (sam_read1 returned %d)", ret);
    }

    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    // fill in how many records we wrote
    if (bc_idx_write_header(out_fp, idx) > 0) {
        FC_LOG_ERROR("Failed to write header to index.");
    }
    fclose(out_fp);
    FC_LOG_INFO("Written %zu/%d records to index.", idx->n_chunks, i);
    bc_idx_destroy(idx);
}

int main_build(int argc, char *argv[])
{
    fc_log_config log_cfg = {
        .level = FC_LOG_INFO,
        .module = "BAMBUILD",
        .target = NULL,
        .color = FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
#ifdef NOTHREADS
    if (args.threads != 1) {
        FC_LOG_WARN("--threads set to %d, but threading not supported by this build.", args.threads);
    }
#endif

    char *index_fname = generate_index_filename(args.bam, NULL);
    index_build(args.bam, index_fname, args.threads, args.chunk_size);
    free(index_fname);

    clock_t end = clock();
    FC_LOG_INFO("Total CPU time: %fs", (double)(end - begin) / CLOCKS_PER_SEC);
    fc_log_destroy();
    return EXIT_SUCCESS;
}
