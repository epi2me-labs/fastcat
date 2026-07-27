#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/fc_accuracy_qscore_2d.h"
#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct {
    fc_accuracy_qscore_2d_sink_opts cfg;
    char *output_root;
} accuracy_qscore_2d_sink;

static FILE *aq2d_open_output(const accuracy_qscore_2d_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (!dir) {
        FC_LOG_FATAL("Failed to create output directory for accuracy_qscore_2d (barcode %zu).", barcode);
    }
    char *path = fc_join_path(dir, "accuracy_qscore_2d.tsv");
    free(dir);
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open accuracy_qscore_2d output file '%s'.", path);
    }
    free(path);
    return fp;
}

static void aq2d_write_matrix(const fc_accuracy_qscore_2d_matrix *matrix, FILE *fp)
{
    fprintf(fp, "concordance_qscore_bin_start\tconcordance_qscore_bin_end\tqscore_bin_start\tqscore_bin_end\tcount\n");
    double width = (matrix->x_upper - matrix->x_lower) / (double)matrix->x_bins;
    for (size_t a = 0; a < matrix->x_bins; a++) {
        double start = matrix->x_lower + (double)a * width;
        double end = (a + 1 == matrix->x_bins) ? matrix->x_upper : start + width;
        for (size_t q = 0; q < matrix->qscore_bins; q++) {
            double count = matrix->counts[q * matrix->x_bins + a];
            if (!isfinite(count) || count < 0.0) {
                count = 0.0;
            }
            double qstart = (double)q * matrix->qscore_width;
            fprintf(fp, "%.10g\t%.10g\t%.1f\t%.1f\t%.12g\n",
                start, end, qstart, qstart + matrix->qscore_width, count);
        }
    }
}

static void aq2d_write_barcode_output(const accuracy_qscore_2d_sink *sink, size_t barcode)
{
    fc_accuracy_qscore_2d_matrix matrix = { 0 };
    if (fc_accuracy_qscore_2d_materialize(sink->cfg.model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize accuracy_qscore_2d matrix for barcode %zu.", barcode);
    }
    FILE *fp = aq2d_open_output(sink, barcode);
    aq2d_write_matrix(&matrix, fp);
    fclose(fp);
    fc_accuracy_qscore_2d_matrix_destroy(&matrix);
}

static void accuracy_qscore_2d_sink_finish(void *ctx)
{
    accuracy_qscore_2d_sink *sink = ctx;
    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (fc_accuracy_qscore_2d_seen(sink->cfg.model, b)) {
                aq2d_write_barcode_output(sink, b);
            }
        }
        return;
    }
    if (fc_accuracy_qscore_2d_seen(sink->cfg.model, 0)) {
        aq2d_write_barcode_output(sink, 0);
    }
}

static void accuracy_qscore_2d_sink_destroy(void *ctx)
{
    accuracy_qscore_2d_sink *sink = ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_accuracy_qscore_2d_sink_create(const fc_accuracy_qscore_2d_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root || !opts->model) {
        return sink;
    }
    accuracy_qscore_2d_sink *state = fc_alloc(1, sizeof(*state), "accuracy_qscore_2d sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);
    sink.ctx = state;
    sink.finish = accuracy_qscore_2d_sink_finish;
    sink.destroy = accuracy_qscore_2d_sink_destroy;
    return sink;
}
