#include "fc_sinks.h"

#include <stdlib.h>
#include <stdio.h>
#include "sink/common.h"

// Coverage is a thin adapter: cov_writer owns most lifecycle and file I/O concerns that would
// otherwise live in the sink/dispatcher layer. Future cleanup could lift output setup/teardown
// (paths, headers, file rotation) into a sink and keep cov_writer focused on per-read updates.
typedef struct {
    cov_writer coverage;
} coverage_sink;

static void coverage_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr) {
    coverage_sink *sink = (coverage_sink *)ctx;
    (void)fname;
    (void)hdr;
    (void)read_ctx;
    if (sink->coverage) {
        coverage_process(sink->coverage, read->bam);
    }
}

static void coverage_sink_destroy(void *ctx) {
    coverage_sink *sink = (coverage_sink *)ctx;
    if (sink->coverage) {
        destroy_coverage_writer(sink->coverage);
    }
    free(sink);
}

fc_sink fc_coverage_sink_create(const fc_coverage_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts || !opts->output_root) return sink;

    char *out_dir = fc_sink_join_path(opts->output_root, "coverage");
    if (!out_dir) return sink;

    cov_writer coverage = init_coverage_writer(
        out_dir,
        opts->per_base,
        opts->use_cigar,
        opts->exclude_flags,
        opts->include_flags,
        opts->by_strand,
        opts->hdr,
        opts->pool,
        opts->bed_files,
        opts->bed_names,
        opts->n_beds,
        opts->thresholds,
        opts->n_thresholds,
        opts->segments,
        opts->n_segments);
    free(out_dir);

    if (!coverage) return sink;
    coverage_sink *state = calloc(1, sizeof(*state));
    if (!state) {
        destroy_coverage_writer(coverage);
        exit(EXIT_FAILURE);
    }
    state->coverage = coverage;
    sink.ctx = state;
    sink.add_read = coverage_sink_add_read;
    sink.destroy = coverage_sink_destroy;
    return sink;
}
