#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "core/fc_stats.h"
#include "sink/common.h"

typedef struct {
    fc_hist_sink_opts cfg; // cached configuration
    char *hist_dir;
    char *output_dir;
} hist_sink;

static void hist_write_stats(const char *hist_dir, const char *plex_dir, size_t barcode, read_stats *stats, const char *type)
{
    char *filepath;
    const char *suff = "hist";

    if (plex_dir == NULL) {
        filepath = fc_asprintf("histogram output path", "%s/%s.%s", hist_dir, type, suff);
    } else {
        char stem_buf[32] = { 0 };
        const char *stem = fc_sink_barcode_stem(barcode, stem_buf, sizeof(stem_buf));
        if (!stem) {
            FC_LOG_FATAL("Could not format barcode stem for histogram output.");
        }
        char *dir = fc_sink_ensure_barcode_dir(plex_dir, barcode);
        if (!dir) {
            FC_LOG_FATAL("Failed to create barcode directory for histograms (barcode %zu).", barcode);
        }
        filepath = fc_asprintf("histogram output path", "%s/%s.%s.%s", dir, stem, type, suff);
        free(dir);
    }

    ensure_parent_dir_exists(filepath);
    FILE *fp = fopen(filepath, "w");
    if (fp == NULL) {
        FC_LOG_FATAL("Could not open file %s for writing", filepath);
    }
    print_stats(stats, false, true, fp);
    fclose(fp);
    free(filepath);
}

static void hist_write_engine(hist_sink *sink, size_t barcode, fc_stats_engine *engine)
{

    if (sink->cfg.unmapped_override) {
        if (engine->length_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->length_stats_unmapped, "length");
        }
        if (engine->qual_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->qual_stats_unmapped, "quality");
        }
    } else {
        if (engine->length_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->length_stats, "length");
        }
        if (engine->qual_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->qual_stats, "quality");
        }
    }

    if (sink->cfg.write_alignment_metrics) {
        if (engine->acc_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->acc_stats, "accuracy");
        }
        if (engine->cov_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->cov_stats, "coverage");
        }
    }

    if (sink->cfg.write_polya && engine->polya_stats) {
        hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->polya_stats, "polya");
    }

    if (!sink->cfg.unmapped_override) {
        const char *len_name = "length.unmap";
        const char *qual_name = "quality.unmap";
        if (engine->length_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->length_stats_unmapped, len_name);
        }
        if (engine->qual_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->qual_stats_unmapped, qual_name);
        }
    }
}

static void hist_sink_finish(void *ctx)
{
    hist_sink *sink = (hist_sink *)ctx;
    for (size_t i = 0; i < sink->cfg.engine_count; ++i) {
        if (sink->cfg.engines[i] != NULL) {
            hist_write_engine(sink, i, sink->cfg.engines[i]);
        }
    }
}

static void hist_sink_destroy(void *ctx)
{
    hist_sink *sink = (hist_sink *)ctx;
    free(sink->hist_dir);
    free(sink->output_dir);
    free(sink);
}

fc_sink fc_hist_sink_create(const fc_hist_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root || !opts->engines || opts->engine_count == 0) {
        return sink;
    }

    hist_sink *state = fc_alloc(1, sizeof(*state), "hist sink");

    state->cfg = *opts;

    if (state->cfg.use_barcodes) {
        // Demultiplexed: write under <output_root>/<barcode>/...; no top-level hist_dir
        state->output_dir = fc_strip_path(opts->output_root);
        state->hist_dir = NULL;
    } else {
        // Single set of histograms under <output_root>/histograms
        char *default_hist = fc_join_path(opts->output_root, "histograms");
        state->hist_dir = fc_strip_path(default_hist);
        free(default_hist);
        state->output_dir = NULL;

        if (state->hist_dir != NULL) {
            int rtn = mkdir_hier(state->hist_dir, true);
            if (rtn == -1) {
                FC_LOG_FATAL("Failed to create histogram directory '%s'.", state->hist_dir);
            }
        }
    }

    sink.ctx = state;
    sink.finish = hist_sink_finish;
    sink.destroy = hist_sink_destroy;
    return sink;
}
