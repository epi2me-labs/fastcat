#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct {
    fc_flagstats_sink_opts cfg; // cached configuration
    char *sample_dup;
    char *region_label_dup;
    FILE *global_fp;
    FILE **bc_fp;
    char *output_root;
} flagstats_sink;

static void flagstats_write_header(FILE *fh, const char *sample)
{
    if (!fh) {
        return;
    }
    fprintf(fh, "ref%s", sample ? "\tsample_name" : "");
    for (size_t i = 0; i < FC_FLAGSTAT_COUNT; ++i) {
        fprintf(fh, "\t%s", fc_flagstat_column_names[i]);
    }
    fputc('\n', fh);
}

static void flagstats_write_line(FILE *fh, const char *sample, const char *label, const size_t *stats)
{
    if (!fh) {
        return;
    }
    fprintf(fh, "%s", label);
    if (sample) {
        fprintf(fh, "\t%s", sample);
    }
    for (size_t i = 0; i < FC_FLAGSTAT_COUNT; ++i) {
        fprintf(fh, "\t%zu", stats[i]);
    }
    fputc('\n', fh);
}

static void flagstats_write_all(flagstats_sink *sink, const fc_alignment_counts *stats, FILE *fp)
{
    if (!stats || !fp) {
        return;
    }

    if (sink->cfg.region_mode) {
        const char *label = sink->cfg.region_label ? sink->cfg.region_label : "region";
        flagstats_write_line(fp, sink->cfg.sample, label, stats->reference_counts[0]);
    } else {
        for (size_t i = 0; i < stats->n_refs; ++i) {
            const char *label = NULL;
            if (sink->cfg.hdr && i < (size_t)sink->cfg.hdr->n_targets) {
                label = sam_hdr_tid2name((sam_hdr_t *)sink->cfg.hdr, (int)i);
            }
            char tmp[32];
            if (!label) {
                FC_SNPRINTF_CHECK(tmp, sizeof(tmp), "%zu", i);
                label = tmp;
            }
            flagstats_write_line(fp, sink->cfg.sample, label, stats->reference_counts[i]);
        }
    }
    if (stats->unmapped_counts) {
        flagstats_write_line(fp, sink->cfg.sample, "*", stats->unmapped_counts);
    }
}

static FILE *flagstats_open_file(flagstats_sink *sink, const char *path)
{
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open file '%s' for writing.", path);
    }
    flagstats_write_header(fp, sink->cfg.sample);
    return fp;
}

static void flagstats_sink_finish(void *ctx)
{
    flagstats_sink *s = (flagstats_sink *)ctx;
    if (!s) {
        return;
    }
    const fc_alignment_counts *global = fc_alignment_summary_global(
        s->cfg.alignment_summary);
    if (s->global_fp && global) {
        flagstats_write_all(s, global, s->global_fp);
    }
    if (s->cfg.use_barcodes) {
        const char *leaf = "flagstats.tsv";
        for (size_t b = 0; b < s->cfg.barcode_count; ++b) {
            const fc_alignment_counts *counts = fc_alignment_summary_barcode(
                s->cfg.alignment_summary,
                b);
            if (!counts) {
                continue;
            }
            if (!s->bc_fp[b]) {
                char *dir = fc_sink_ensure_barcode_dir(s->output_root, b);
                if (!dir) {
                    FC_LOG_FATAL("Failed to create barcode directory for flagstats (barcode %zu).", b);
                }
                char *path = fc_join_path(dir, leaf);
                free(dir);
                s->bc_fp[b] = flagstats_open_file(s, path);
                free(path);
            }
            flagstats_write_all(s, counts, s->bc_fp[b]);
        }
    }
}

static void flagstats_sink_destroy(void *ctx)
{
    flagstats_sink *sink = (flagstats_sink *)ctx;
    if (sink->global_fp) {
        fclose(sink->global_fp);
    }
    if (sink->bc_fp) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->bc_fp[i]) {
                fclose(sink->bc_fp[i]);
            }
        }
    }
    free(sink->sample_dup);
    free(sink->region_label_dup);
    free(sink->bc_fp);
    free(sink->output_root);
    free(sink);
}

fc_sink fc_flagstats_sink_create(const fc_flagstats_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root || !opts->alignment_summary) {
        return sink;
    }

    flagstats_sink *state = fc_alloc(1, sizeof(*state), "flagstats sink");
    state->cfg = *opts;
    state->cfg.sample = NULL;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->bc_fp = fc_alloc(state->cfg.barcode_count, sizeof(FILE *), "flagstats barcode files");
    state->output_root = fc_strip_path(opts->output_root);
    if (opts->sample) {
        state->sample_dup = fc_strdup(opts->sample, "sample name for flagstats sink");
        state->cfg.sample = state->sample_dup;
    }
    if (opts->region_label) {
        state->region_label_dup = fc_strdup(opts->region_label, "region label for flagstats sink");
        state->cfg.region_label = state->region_label_dup;
    }
    if (state->output_root) {
        char *path = fc_join_path(state->output_root, "flagstats.tsv");
        state->global_fp = flagstats_open_file(state, path);
        free(path);
    }

    sink.ctx = state;
    sink.finish = flagstats_sink_finish;
    sink.destroy = flagstats_sink_destroy;
    return sink;
}
