#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../fc_log.h"
#include "common.h"
#include "fc_sinks.h"
#include "kh_counter.h"
#include "sink/common.h"

typedef struct {
    fc_counter_sink_opts cfg; // cached configuration
    FILE *runids_global;
    FILE *basecallers_global;
    FILE **runids_bc;
    FILE **basecallers_bc;
    char *sample; // duplicated if provided
    char *output_root; // stripped output root
    kh_counter_t *run_ids_global_ctr;
    kh_counter_t *basecallers_global_ctr;
    kh_counter_t **run_ids_bc_ctr;
    kh_counter_t **basecallers_bc_ctr;
} counter_sink;

static void destroy_counter(kh_counter_t **ctr)
{
    if (*ctr) {
        kh_counter_destroy(*ctr);
        *ctr = NULL;
    }
}

static void counter_reset(counter_sink *sink)
{
    destroy_counter(&sink->run_ids_global_ctr);
    destroy_counter(&sink->basecallers_global_ctr);
    if (sink->cfg.use_barcodes) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            destroy_counter(&sink->run_ids_bc_ctr[i]);
            destroy_counter(&sink->basecallers_bc_ctr[i]);
        }
    }
}

static kh_counter_t *ensure_counter(kh_counter_t **arr, size_t idx)
{
    if (!arr[idx]) {
        arr[idx] = kh_counter_init();
    }
    return arr[idx];
}

static FILE *counter_open_file(counter_sink *sink, const char *path, const char *value_label)
{
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open counter file '%s' for writing.", path);
    }
    fprintf(fp, "filename\t");
    if (sink->sample) {
        fprintf(fp, "sample_name\t");
    }
    fprintf(fp, "%s\tcount\n", value_label);
    return fp;
}

static FILE *ensure_bc_file(counter_sink *sink, FILE **files, const char *leaf, const char *value_label, size_t barcode)
{
    if (!sink->cfg.use_barcodes || barcode >= sink->cfg.barcode_count) {
        return NULL;
    }
    if (files[barcode]) {
        return files[barcode];
    }

    char *dir = fc_sink_ensure_barcode_dir(sink->output_root, barcode);
    if (!dir) {
        FC_LOG_ERROR("Failed to create barcode directory for %s (barcode %zu).", leaf, barcode);
        return NULL;
    }
    char *path = fc_join_path(dir, leaf);
    free(dir);
    FILE *fp = counter_open_file(sink, path, value_label);
    free(path);
    files[barcode] = fp;
    return fp;
}

static void counter_sink_begin_file(void *ctx, const char *fname)
{
    counter_sink *sink = (counter_sink *)ctx;
    (void)fname;
    counter_reset(sink);
}

static void counter_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
    const char *fname, const sam_hdr_t *hdr)
{
    counter_sink *sink = (counter_sink *)ctx;
    (void)fname;
    (void)hdr;
    size_t barcode = (sink->cfg.use_barcodes && read_ctx) ? read_ctx->barcode : 0;

    const char *runid = read->meta->runid ? read->meta->runid : "";
    const char *basecaller = read->meta->basecaller ? read->meta->basecaller : "";
    kh_counter_increment(ensure_counter(&sink->run_ids_global_ctr, 0), (char *)runid);
    kh_counter_increment(ensure_counter(&sink->basecallers_global_ctr, 0), (char *)basecaller);

    if (sink->cfg.use_barcodes && barcode < sink->cfg.barcode_count) {
        kh_counter_increment(ensure_counter(sink->run_ids_bc_ctr, barcode), (char *)runid);
        kh_counter_increment(ensure_counter(sink->basecallers_bc_ctr, barcode), (char *)basecaller);
    }
}

static void counter_sink_end_file(void *ctx, const char *fname, const fc_file_stats *stats)
{
    counter_sink *sink = (counter_sink *)ctx;
    (void)stats;
    if (sink->runids_global && sink->run_ids_global_ctr) {
        for (khiter_t k = 0; k < kh_end(sink->run_ids_global_ctr); ++k) {
            if (kh_exist(sink->run_ids_global_ctr, k)) {
                fprintf(sink->runids_global, "%s\t", fname);
                if (sink->sample) {
                    fprintf(sink->runids_global, "%s\t", sink->sample);
                }
                fprintf(sink->runids_global, "%s\t%d\n", kh_key(sink->run_ids_global_ctr, k), kh_val(sink->run_ids_global_ctr, k));
            }
        }
    }
    if (sink->basecallers_global && sink->basecallers_global_ctr) {
        for (khiter_t k = 0; k < kh_end(sink->basecallers_global_ctr); ++k) {
            if (kh_exist(sink->basecallers_global_ctr, k)) {
                fprintf(sink->basecallers_global, "%s\t", fname);
                if (sink->sample) {
                    fprintf(sink->basecallers_global, "%s\t", sink->sample);
                }
                fprintf(sink->basecallers_global, "%s\t%d\n", kh_key(sink->basecallers_global_ctr, k), kh_val(sink->basecallers_global_ctr, k));
            }
        }
    }
    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; ++b) {
            if (sink->run_ids_bc_ctr[b]) {
                FILE *fp = ensure_bc_file(sink, sink->runids_bc, "runids.tsv", "run_id", b);
                for (khiter_t k = 0; k < kh_end(sink->run_ids_bc_ctr[b]); ++k) {
                    if (kh_exist(sink->run_ids_bc_ctr[b], k)) {
                        fprintf(fp, "%s\t", fname);
                        if (sink->sample) {
                            fprintf(fp, "%s\t", sink->sample);
                        }
                        fprintf(fp, "%s\t%d\n", kh_key(sink->run_ids_bc_ctr[b], k), kh_val(sink->run_ids_bc_ctr[b], k));
                    }
                }
            }
            if (sink->basecallers_bc_ctr[b]) {
                FILE *fp = ensure_bc_file(sink, sink->basecallers_bc, "basecallers.tsv", "basecaller", b);
                for (khiter_t k = 0; k < kh_end(sink->basecallers_bc_ctr[b]); ++k) {
                    if (kh_exist(sink->basecallers_bc_ctr[b], k)) {
                        fprintf(fp, "%s\t", fname);
                        if (sink->sample) {
                            fprintf(fp, "%s\t", sink->sample);
                        }
                        fprintf(fp, "%s\t%d\n", kh_key(sink->basecallers_bc_ctr[b], k), kh_val(sink->basecallers_bc_ctr[b], k));
                    }
                }
            }
        }
    }
}

static void counter_sink_destroy(void *ctx)
{
    counter_sink *sink = (counter_sink *)ctx;
    if (sink->runids_global) {
        fclose(sink->runids_global);
    }
    if (sink->basecallers_global) {
        fclose(sink->basecallers_global);
    }
    if (sink->runids_bc) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->runids_bc[i]) {
                fclose(sink->runids_bc[i]);
            }
        }
    }
    if (sink->basecallers_bc) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->basecallers_bc[i]) {
                fclose(sink->basecallers_bc[i]);
            }
        }
    }
    destroy_counter(&sink->run_ids_global_ctr);
    destroy_counter(&sink->basecallers_global_ctr);
    if (sink->run_ids_bc_ctr) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            destroy_counter(&sink->run_ids_bc_ctr[i]);
        }
    }
    if (sink->basecallers_bc_ctr) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            destroy_counter(&sink->basecallers_bc_ctr[i]);
        }
    }
    free(sink->runids_bc);
    free(sink->basecallers_bc);
    free(sink->run_ids_bc_ctr);
    free(sink->basecallers_bc_ctr);
    free(sink->sample);
    free(sink->output_root);
    free(sink);
}

fc_sink fc_counter_sink_create(const fc_counter_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root) {
        return sink;
    }

    counter_sink *state = fc_alloc(1, sizeof(*state), "counter sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);
    if (opts->sample) {
        state->sample = fc_strdup(opts->sample, "sample name for counter sink");
        state->cfg.sample = state->sample;
    }
    char *runids_path = fc_join_path(opts->output_root, "runids.tsv");
    state->runids_global = counter_open_file(state, runids_path, "run_id");
    free(runids_path);

    char *basecallers_path = fc_join_path(opts->output_root, "basecallers.tsv");
    state->basecallers_global = counter_open_file(state, basecallers_path, "basecaller");
    free(basecallers_path);

    state->run_ids_bc_ctr = fc_alloc(state->cfg.barcode_count, sizeof(kh_counter_t *), "barcode run-id counters");
    state->basecallers_bc_ctr = fc_alloc(state->cfg.barcode_count, sizeof(kh_counter_t *), "barcode basecaller counters");
    state->runids_bc = fc_alloc(state->cfg.barcode_count, sizeof(FILE *), "barcode run-id files");
    state->basecallers_bc = fc_alloc(state->cfg.barcode_count, sizeof(FILE *), "barcode basecaller files");
    state->run_ids_global_ctr = NULL;
    state->basecallers_global_ctr = NULL;

    sink.ctx = state;
    sink.begin_file = counter_sink_begin_file;
    sink.add_read = counter_sink_add_read;
    sink.end_file = counter_sink_end_file;
    sink.destroy = counter_sink_destroy;
    return sink;
}
