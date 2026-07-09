#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_run_summary.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct {
    fc_run_summary_sink_opts cfg;
    char *output_root;
} run_summary_sink;

static double run_summary_stats_total_count(const read_stats *stats)
{
    if (!stats) {
        return 0.0;
    }
    double total = 0.0;
    for (size_t i = 0; i < stats->n; ++i) {
        total += (double)stats->counts[i];
    }
    return total;
}

static const read_stats *run_summary_length_stats(const fc_stats_engine *engine, bool unmapped_override)
{
    if (!engine) {
        return NULL;
    }
    if (unmapped_override) {
        return engine->length_stats_unmapped;
    }
    if (engine->length_stats && run_summary_stats_total_count(engine->length_stats) > 0.0) {
        return engine->length_stats;
    }
    if (engine->length_stats_unmapped && run_summary_stats_total_count(engine->length_stats_unmapped) > 0.0) {
        return engine->length_stats_unmapped;
    }
    return engine->length_stats;
}

static FILE *run_summary_open_output(const run_summary_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (NULL == dir) {
        FC_LOG_FATAL("Failed to create output directory for run_summary (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "run_summary.tsv");
    free(dir);
    ensure_parent_dir_exists(path);

    FILE *fp = fopen(path, "w");
    if (NULL == fp) {
        FC_LOG_FATAL("Cannot open run_summary output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

static void run_summary_write_row(FILE *fp, const fc_run_summary_row *row)
{
    fprintf(fp,
        "sample_id\tinput_type\tbarcode_id\tbarcode_alias\trun_id\tbasecaller_model\tflow_cell_id\trun_date\t"
        "run_duration_s\tread_count\tbase_count\tn50_bp\tl50\n");
    fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
        row->sample_id ? row->sample_id : "",
        row->input_type ? row->input_type : "",
        row->barcode_id ? row->barcode_id : "",
        row->barcode_alias ? row->barcode_alias : "",
        row->run_id ? row->run_id : "",
        row->basecaller_model ? row->basecaller_model : "",
        row->flow_cell_id ? row->flow_cell_id : "",
        row->run_date ? row->run_date : "");
    if (isfinite(row->run_duration_s)) {
        fprintf(fp, "%.12g", row->run_duration_s);
    } else {
        fputs("nan", fp);
    }
    fprintf(fp, "\t%zu\t%zu\t%zu\t%zu\n",
        row->read_count,
        row->base_count,
        row->n50_bp,
        row->l50);
}

static void run_summary_write_output(const run_summary_sink *sink, size_t barcode)
{
    if (!fc_run_summary_seen(sink->cfg.model, barcode)) {
        return;
    }

    char barcode_buf[32] = { 0 };
    const char *barcode_id = sink->cfg.use_barcodes
        ? fc_sink_barcode_stem(barcode, barcode_buf, sizeof(barcode_buf))
        : "all";
    const read_stats *length_stats = NULL;
    if (sink->cfg.engines && barcode < sink->cfg.engine_count) {
        length_stats = run_summary_length_stats(sink->cfg.engines[barcode], sink->cfg.unmapped_override);
    }

    fc_run_summary_row row = { 0 };
    if (fc_run_summary_materialize(
            sink->cfg.model,
            barcode,
            (sink->cfg.engines && barcode < sink->cfg.engine_count) ? sink->cfg.engines[barcode] : NULL,
            sink->cfg.bam_input,
            length_stats,
            sink->cfg.sample,
            barcode_id,
            &row)
        != 0) {
        FC_LOG_FATAL("Failed to materialize run_summary row for barcode %zu.", barcode);
    }

    FILE *fp = run_summary_open_output(sink, barcode);
    run_summary_write_row(fp, &row);
    fclose(fp);
}

static void run_summary_sink_finish(void *ctx)
{
    run_summary_sink *sink = (run_summary_sink *)ctx;
    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; ++b) {
            run_summary_write_output(sink, b);
        }
    } else {
        run_summary_write_output(sink, 0);
    }
}

static void run_summary_sink_destroy(void *ctx)
{
    run_summary_sink *sink = (run_summary_sink *)ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_run_summary_sink_create(const fc_run_summary_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (NULL == opts || NULL == opts->output_root || NULL == opts->model) {
        return sink;
    }

    run_summary_sink *state = fc_alloc(1, sizeof(*state), "run_summary sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);

    sink.ctx = state;
    sink.finish = run_summary_sink_finish;
    sink.destroy = run_summary_sink_destroy;
    return sink;
}
