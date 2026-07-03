#include <stdio.h>
#include <stdlib.h>

#include "core/fc_channel_summary.h"
#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct {
    fc_channel_summary_sink_opts cfg;
    char *output_root;
} channel_sink;

static FILE *channel_open_output(const channel_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (NULL == dir) {
        FC_LOG_FATAL("Failed to create output directory for channel_summary (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "channel_summary.tsv");
    free(dir);
    ensure_parent_dir_exists(path);

    FILE *fp = fopen(path, "w");
    if (NULL == fp) {
        FC_LOG_FATAL("Cannot open channel_summary output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

static void channel_write_table(const fc_channel_summary_table *table, FILE *fp)
{
    fprintf(fp, "channel\tmux\tread_count\tbase_count\n");
    for (size_t i = 0; i < table->n; i++) {
        const fc_channel_summary_row *row = &table->rows[i];
        fprintf(fp, "%zu\t%zu\t%zu\t%zu\n",
            row->channel,
            row->mux,
            row->read_count,
            row->base_count);
    }
}

static void channel_write_barcode_output(const channel_sink *sink, size_t barcode)
{
    fc_channel_summary_table table = { 0 };
    if (fc_channel_summary_materialize_sparse(sink->cfg.model, barcode, &table) != 0) {
        FC_LOG_FATAL("Failed to materialize channel_summary table for barcode %zu.", barcode);
    }

    FILE *fp = channel_open_output(sink, barcode);
    channel_write_table(&table, fp);
    fclose(fp);
    fc_channel_summary_table_destroy(&table);
}

static void channel_sink_finish(void *ctx)
{
    channel_sink *sink = (channel_sink *)ctx;
    fc_channel_summary_skip_counts skipped = fc_channel_summary_skipped_counts(sink->cfg.model);
    if (skipped.missing_channel || skipped.invalid_channel || skipped.invalid_mux) {
        FC_LOG_WARN(
            "Skipped %zu reads without channel, %zu with invalid channel and %zu with invalid mux for channel_summary.",
            skipped.missing_channel,
            skipped.invalid_channel,
            skipped.invalid_mux);
    }

    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (fc_channel_summary_seen(sink->cfg.model, b)) {
                channel_write_barcode_output(sink, b);
            }
        }
    } else {
        channel_write_barcode_output(sink, 0);
    }
}

static void channel_sink_destroy(void *ctx)
{
    channel_sink *sink = (channel_sink *)ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_channel_summary_sink_create(const fc_channel_summary_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (NULL == opts || NULL == opts->output_root || NULL == opts->model) {
        return sink;
    }

    channel_sink *state = fc_alloc(1, sizeof(*state), "channel sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);

    sink.ctx = state;
    sink.finish = channel_sink_finish;
    sink.destroy = channel_sink_destroy;
    return sink;
}
