#include <stdio.h>
#include <stdlib.h>

#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

// Coverage is a thin adapter: cov_writer owns most lifecycle and file I/O concerns that would
// otherwise live in the sink/dispatcher layer. Future cleanup could lift output setup/teardown
// (paths, headers, file rotation) into a sink and keep cov_writer focused on per-read updates.
typedef struct {
    fc_coverage_sink_opts cfg; // cached configuration
    cov_writer global_writer;
    cov_writer *barcode_writers;
    char *output_root; // stripped output root
} coverage_sink;

static cov_writer coverage_create_writer(coverage_sink *sink, size_t barcode)
{
    char *out_dir = NULL;
    if (sink->cfg.write_traces || sink->cfg.write_tables) {
        char *dir = fc_sink_ensure_barcode_dir(sink->output_root, barcode);
        if (!dir) {
            FC_LOG_FATAL("Failed to create barcode directory under '%s' (barcode %zu).",
                sink->output_root, barcode);
        }
        out_dir = fc_join_path(dir, "coverage");
        free(dir);
    }
    cov_writer w = init_coverage_writer(
        out_dir,
        sink->cfg.write_traces,
        sink->cfg.write_tables,
        sink->cfg.use_cigar,
        sink->cfg.exclude_flags,
        sink->cfg.include_flags,
        sink->cfg.by_strand,
        sink->cfg.hdr,
        sink->cfg.pool,
        sink->cfg.bed_files,
        sink->cfg.bed_names,
        sink->cfg.n_beds,
        sink->cfg.thresholds,
        sink->cfg.n_thresholds,
        sink->cfg.segments,
        sink->cfg.n_segments,
        fc_coverage_summaries_prepare(sink->cfg.summaries, barcode));
    free(out_dir);
    return w;
}

static void coverage_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
    const char *fname, const sam_hdr_t *hdr)
{
    coverage_sink *sink = (coverage_sink *)ctx;
    (void)fname;
    (void)hdr;
    if (sink->global_writer) {
        coverage_process(sink->global_writer, read->bam);
    }
    if (sink->cfg.use_barcodes && read_ctx) {
        size_t b = read_ctx->barcode;
        if (b < sink->cfg.barcode_count) {
            if (!sink->barcode_writers[b]) {
                sink->barcode_writers[b] = coverage_create_writer(sink, b);
            }
            if (sink->barcode_writers[b]) {
                coverage_process(sink->barcode_writers[b], read->bam);
            }
        }
    }
}

static void coverage_sink_destroy(void *ctx)
{
    coverage_sink *sink = (coverage_sink *)ctx;
    if (sink->global_writer) {
        destroy_coverage_writer(sink->global_writer);
    }
    if (sink->barcode_writers) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->barcode_writers[i]) {
                destroy_coverage_writer(sink->barcode_writers[i]);
            }
        }
        free(sink->barcode_writers);
    }
    free(sink->output_root);
    free(sink);
}

static void coverage_sink_finish(void *ctx)
{
    coverage_sink *sink = (coverage_sink *)ctx;
    if (sink->global_writer) {
        destroy_coverage_writer(sink->global_writer);
        sink->global_writer = NULL;
    }
    for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
        if (sink->barcode_writers[i]) {
            destroy_coverage_writer(sink->barcode_writers[i]);
            sink->barcode_writers[i] = NULL;
        }
    }
}

fc_sink fc_coverage_sink_create(const fc_coverage_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || ((opts->write_traces || opts->write_tables) && !opts->output_root)) {
        return sink;
    }

    coverage_sink *state = fc_alloc(1, sizeof(*state), "coverage sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    if (state->cfg.write_traces || state->cfg.write_tables) {
        state->output_root = fc_strip_path(opts->output_root);
    }
    if ((state->cfg.write_traces || state->cfg.write_tables) && !state->output_root) {
        free(state);
        return sink;
    }
    if (!state->cfg.use_barcodes) {
        char *out_dir = NULL;
        if (state->cfg.write_traces || state->cfg.write_tables) {
            out_dir = fc_join_path(state->output_root, "coverage");
        }
        cov_writer coverage = init_coverage_writer(
            out_dir,
            state->cfg.write_traces,
            state->cfg.write_tables,
            state->cfg.use_cigar,
            state->cfg.exclude_flags,
            state->cfg.include_flags,
            state->cfg.by_strand,
            state->cfg.hdr,
            state->cfg.pool,
            state->cfg.bed_files,
            state->cfg.bed_names,
            state->cfg.n_beds,
            state->cfg.thresholds,
            state->cfg.n_thresholds,
            state->cfg.segments,
            state->cfg.n_segments,
            fc_coverage_summaries_prepare(state->cfg.summaries, 0));
        free(out_dir);
        state->global_writer = coverage;
        if (!state->global_writer) {
            free(state->output_root);
            free(state);
            return sink;
        }
    }
    state->barcode_writers = fc_alloc(state->cfg.barcode_count, sizeof(cov_writer), "coverage barcode writers");
    sink.ctx = state;
    sink.add_read = coverage_sink_add_read;
    sink.finish = coverage_sink_finish;
    sink.destroy = coverage_sink_destroy;
    return sink;
}
