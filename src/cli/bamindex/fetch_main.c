// bamindex fetch program

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "core/fc_argp.h"
#include "core/fc_bam_utils.h"
#include "core/fc_common.h"
#include "core/fc_log.h"
#include "index.h"

static char doc[] = "bamindex fetch -- fetch records from a BAM according to an index.\
\vThe program simply will fetch a batch of records from a BAM file"
                    "using and index and a chunk ID. Output is written as uncompressed "
                    "BAM to stdout.";
static char args_doc[] = "<reads.bam.bci>";
static struct argp_option options[] = {
    { 0, 0, 0, 0,
        "General options:", 0 },
    { "threads", 't', "THREADS", 0,
        "Number of threads for BAM processing.", 0 },
    { "chunk", 'c', "SIZE", 0,
        "Chunk index to retrieve.", 0 },
    { 0 }
};

typedef struct arguments {
    const char *bam;
    char *index;
    int chunk_idx;
    int threads;
} arguments_t;

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
            argp_error(state, "Invalid --chunk value: '%s'.", arg);
        }
        if (int_value < 0) {
            argp_error(state, "--chunk must be >= 0.");
        }
        arguments->chunk_idx = int_value;
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
    args.index = NULL;
    args.threads = 1;
    args.chunk_idx = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    args.index = generate_index_filename(args.bam, args.index);
    return args;
}

void index_fetch(const char *bam_fname, const char *index_fname, int chunk, int threads)
{
    htsFile *fp = hts_open(bam_fname, "r");
    if (fp == NULL) {
        FC_LOG_FATAL("Could not open %s", bam_fname);
    }
    fc_require_bam_input(fp, bam_fname, "bamindex fetch");
    bam_hdr_t *h = sam_hdr_read(fp);
    if (h == NULL) {
        FC_LOG_FATAL("Could not read BAM header from %s", bam_fname);
    }
    htsThreadPool p = { NULL, 0 };
    if (threads > 1) {
        p.pool = hts_tpool_init(threads);
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
    }

    struct stat st;
    if (stat(index_fname, &st) != 0) {
        fprintf(stderr, "Cannot open index file %s\n", index_fname);
        exit(EXIT_FAILURE);
    }
    FILE *idx_fp = fopen(index_fname, "rb");
    if (!idx_fp) {
        FC_LOG_FATAL("Cannot open index file: %s.", index_fname);
    }
    bc_idx_t *idx;
    if ((idx = bc_idx_read(idx_fp)) == NULL) {
        fclose(idx_fp);
        FC_LOG_FATAL("Couldn't read index file: %s.", index_fname);
    }
    fclose(idx_fp);
    if ((size_t)chunk >= idx->n_chunks) {
        FC_LOG_ERROR("Requested chunk %d is out of range for index with %zu chunk(s).", chunk, idx->n_chunks);
        bc_idx_destroy(idx);
        bam_hdr_destroy(h);
        hts_close(fp);
        if (p.pool) {
            hts_tpool_destroy(p.pool);
        }
        exit(EXIT_FAILURE);
    }
    bc_rec_t rec = idx->recs[chunk];

    FC_LOG_INFO("Starting from: %zu %s", rec.file_offset, rec.qname);
    FC_LOG_INFO("Reading %zu records from bam.", idx->chunk_size);

    if (bgzf_seek(fp->fp.bgzf, rec.file_offset, SEEK_SET) != 0) {
        FC_LOG_FATAL("Failed to seek to first record.");
    }

    size_t written = 0;
    bam1_t *b = bam_init1();
    htsFile *out_fp;
    if ((out_fp = hts_open("-", "wb0")) == 0) {
        FC_LOG_FATAL("Failed to open standard output for writing.");
    }

    // TODO: fill in the NULLs here
    if (sam_hdr_add_pg(h, "bamindex.fetch", "VN", argp_program_version, NULL, NULL, NULL) != 0) {
        FC_LOG_FATAL("Failed to add PG line to the header.");
    }
    if (sam_hdr_write(out_fp, h) != 0) {
        FC_LOG_FATAL("Failed to write the SAM header.");
    }
    while ((sam_read1(fp, h, b) >= 0) && (written < (idx->chunk_size))) {
        if ((sam_write1(out_fp, h, b) < 0)) {
            FC_LOG_FATAL("Failed to write output record.");
        }
        written++;
    }
    hts_close(out_fp);

    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);
    if (p.pool) { // must be after fp
        hts_tpool_destroy(p.pool);
    }

    bc_idx_destroy(idx);
    FC_LOG_INFO("Written %zu records to output.", written);
}

int main_fetch(int argc, char *argv[])
{
    clock_t begin = clock();
    fc_log_config log_cfg = {
        .level = FC_LOG_INFO,
        .module = "BAMFETCH",
        .target = NULL,
        .color = FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    arguments_t args = parse_arguments(argc, argv);
    index_fetch(args.bam, args.index, args.chunk_idx, args.threads);
    free(args.index);
    clock_t end = clock();
    FC_LOG_INFO("Total CPU time: %fs", (double)(end - begin) / CLOCKS_PER_SEC);
    fc_log_destroy();
    return EXIT_SUCCESS;
}
