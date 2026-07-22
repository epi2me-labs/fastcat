#include <stdio.h>
#include <stdlib.h>

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "core/fc_time_bin.h"
#include "sink/common.h"

typedef struct {
    fc_time_bin_sink_opts cfg;
    char *output_root;
} time_bin_sink;

static FILE *time_bin_open_output(const time_bin_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (NULL == dir) {
        FC_LOG_FATAL("Failed to create output directory for time_bin_summary (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "time_bin_summary.tsv");
    free(dir);
    ensure_parent_dir_exists(path);

    FILE *fp = fopen(path, "w");
    if (NULL == fp) {
        FC_LOG_FATAL("Cannot open time_bin_summary output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

static void time_bin_write_table(const fc_time_bin_table *table, FILE *fp)
{
    fprintf(fp,
        "time_bin_start_s\ttime_bin_end_s\tread_count\tbase_count\t"
        "cumulative_read_count\tcumulative_base_count\t"
        "length_mean\tlength_q05\tlength_q10\tlength_q25\tlength_q50\tlength_q75\tlength_q90\tlength_q95\t"
        "qscore_mean\tqscore_q05\tqscore_q10\tqscore_q25\tqscore_q50\tqscore_q75\tqscore_q90\tqscore_q95\t"
        "translocation_speed_mean\ttranslocation_speed_q05\ttranslocation_speed_q10\ttranslocation_speed_q25\t"
        "translocation_speed_q50\ttranslocation_speed_q75\ttranslocation_speed_q90\ttranslocation_speed_q95\n");

    for (size_t i = 0; i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        fprintf(fp,
            "%.12g\t%.12g\t%zu\t%zu\t%zu\t%zu\t"
            "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t"
            "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t"
            "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
            row->time_bin_start_s,
            row->time_bin_end_s,
            row->read_count,
            row->base_count,
            row->cumulative_read_count,
            row->cumulative_base_count,
            row->length_mean,
            row->length_quantiles.q05,
            row->length_quantiles.q10,
            row->length_quantiles.q25,
            row->length_quantiles.q50,
            row->length_quantiles.q75,
            row->length_quantiles.q90,
            row->length_quantiles.q95,
            row->qscore_mean,
            row->qscore_quantiles.q05,
            row->qscore_quantiles.q10,
            row->qscore_quantiles.q25,
            row->qscore_quantiles.q50,
            row->qscore_quantiles.q75,
            row->qscore_quantiles.q90,
            row->qscore_quantiles.q95,
            row->speed_mean,
            row->speed_quantiles.q05,
            row->speed_quantiles.q10,
            row->speed_quantiles.q25,
            row->speed_quantiles.q50,
            row->speed_quantiles.q75,
            row->speed_quantiles.q90,
            row->speed_quantiles.q95);
    }
}

static void time_bin_write_barcode_output(const time_bin_sink *sink, size_t barcode)
{
    fc_time_bin_table table = { 0 };
    if (fc_time_bin_materialize(sink->cfg.model, barcode, &table) != 0) {
        FC_LOG_FATAL("Failed to materialize time_bin_summary table for barcode %zu.", barcode);
    }

    FILE *fp = time_bin_open_output(sink, barcode);
    time_bin_write_table(&table, fp);
    fclose(fp);
    fc_time_bin_table_destroy(&table);
}

static void time_bin_sink_finish(void *ctx)
{
    time_bin_sink *sink = (time_bin_sink *)ctx;
    fc_time_bin_skip_counts skipped = fc_time_bin_skipped_counts(sink->cfg.model);
    if (skipped.missing || skipped.invalid) {
        FC_LOG_WARN("Skipped %zu reads without start_time and %zu reads with invalid start_time for time_bin_summary.",
            skipped.missing, skipped.invalid);
    }
    if (skipped.duration_missing || skipped.duration_invalid) {
        FC_LOG_WARN("Skipped %zu reads without du and %zu reads with invalid du for time_bin_summary speed.",
            skipped.duration_missing, skipped.duration_invalid);
    }

    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (fc_time_bin_seen(sink->cfg.model, b)) {
                time_bin_write_barcode_output(sink, b);
            }
        }
    } else {
        time_bin_write_barcode_output(sink, 0);
    }
}

static void time_bin_sink_destroy(void *ctx)
{
    time_bin_sink *sink = (time_bin_sink *)ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_time_bin_sink_create(const fc_time_bin_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (NULL == opts || NULL == opts->output_root || NULL == opts->model) {
        return sink;
    }

    time_bin_sink *state = fc_alloc(1, sizeof(*state), "time_bin sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);

    sink.ctx = state;
    sink.finish = time_bin_sink_finish;
    sink.destroy = time_bin_sink_destroy;
    return sink;
}
