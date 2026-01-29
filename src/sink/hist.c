#include "fc_sinks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "common.h"
#include "sink/common.h"
#include "stats.h"

typedef struct {
    char *hist_dir;
    char *output_dir;
    bool use_barcodes;
    bool write_unmapped;
    bool write_alignment_metrics;
    bool write_polya;
    bool unmapped_is_primary;
    size_t engine_count;
    fc_stats_engine **engines;
} hist_sink;

static void hist_write_stats(const char *hist_dir, const char *plex_dir, size_t barcode, read_stats *stats, const char *type) {
    char *filepath;
    const char *suff = "hist";
    if (plex_dir == NULL) {
        filepath = calloc(strlen(hist_dir) + strlen(type) + strlen(suff) + 3, sizeof(char));
        sprintf(filepath, "%s/%s.%s", hist_dir, type, suff);
    } else {
        char *path;
        if (barcode == 0) {
            path = calloc(strlen(plex_dir) + 15, sizeof(char));
            sprintf(path, "%s/unclassified/", plex_dir);
            filepath = calloc(strlen(path) + strlen(type) + strlen(suff) + 15, sizeof(char));
            sprintf(filepath, "%sunclassified.%s.%s", path, type, suff);
        } else {
            path = calloc(strlen(plex_dir) + 14, sizeof(char));
            sprintf(path, "%s/barcode%04lu/", plex_dir, barcode);
            filepath = calloc(strlen(path) + strlen(type) + strlen(suff) + 14, sizeof(char));
            sprintf(filepath, "%sbarcode%04lu.%s.%s", path, barcode, type, suff);
        }
        free(path);
    }
    ensure_parent_dir_exists(filepath);
    FILE *fp = fopen(filepath, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filepath);
        exit(1);
    }
    print_stats(stats, false, true, fp);
    fclose(fp);
    free(filepath);
}

static void hist_write_engine(hist_sink *sink, size_t barcode, fc_stats_engine *engine) {
    if (sink->unmapped_is_primary) {
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

    if (sink->write_alignment_metrics) {
        if (engine->acc_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->acc_stats, "accuracy");
        }
        if (engine->cov_stats) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->cov_stats, "coverage");
        }
    }

    if (sink->write_polya && engine->polya_stats) {
        hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->polya_stats, "polya");
    }

    if (sink->write_unmapped) {
        const char *len_name = sink->unmapped_is_primary ? "length" : "length.unmap";
        const char *qual_name = sink->unmapped_is_primary ? "quality" : "quality.unmap";
        if (engine->length_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->length_stats_unmapped, len_name);
        }
        if (engine->qual_stats_unmapped) {
            hist_write_stats(sink->hist_dir, sink->output_dir, barcode, engine->qual_stats_unmapped, qual_name);
        }
    }
}

static void hist_sink_finish(void *ctx) {
    hist_sink *sink = (hist_sink *)ctx;
    for (size_t i = 0; i < sink->engine_count; ++i) {
        if (sink->engines[i] != NULL) {
            hist_write_engine(sink, i, sink->engines[i]);
        }
    }
}

static void hist_sink_destroy(void *ctx) {
    hist_sink *sink = (hist_sink *)ctx;
    (void)sink->engines;
    free(sink->hist_dir);
    free(sink->output_dir);
    free(sink);
}

fc_sink fc_hist_sink_create(const fc_hist_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts || !opts->output_root || !opts->engines || opts->engine_count == 0) return sink;

    hist_sink *state = calloc(1, sizeof(*state));

    if (opts->use_barcodes) {
        // Demultiplexed: write under <output_root>/demux/...; no top-level hist_dir
        char *default_demux = fc_sink_join_path(opts->output_root, "demux");
        state->output_dir = fc_sink_strip_path(default_demux);
        free(default_demux);
        state->hist_dir = NULL;
        state->use_barcodes = true;
    } else {
        // Single set of histograms under <output_root>/histograms
        char *default_hist = fc_sink_join_path(opts->output_root, "histograms");
        state->hist_dir = fc_sink_strip_path(default_hist);
        free(default_hist);
        state->output_dir = NULL;
        state->use_barcodes = false;

        if (state->hist_dir != NULL) {
            int rtn = mkdir_hier(state->hist_dir);
            if (rtn == -1 && errno != EEXIST) {
                fprintf(stderr,
                    "Error: Cannot create output directory '%s'. Check location is writeable and directory does not exist.\n",
                    state->hist_dir);
                exit(EXIT_FAILURE);
            }
        }
    }

    state->write_unmapped = opts->write_unmapped;
    state->write_alignment_metrics = opts->write_alignment_metrics;
    state->write_polya = opts->write_polya;
    state->unmapped_is_primary = opts->unmapped_is_primary;
    state->engine_count = opts->engine_count;
    state->engines = opts->engines;

    sink.ctx = state;
    sink.finish = hist_sink_finish;
    sink.destroy = hist_sink_destroy;
    return sink;
}
