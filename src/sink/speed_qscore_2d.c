#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "core/fc_speed_qscore_2d.h"
#include "sink/common.h"

typedef struct {
    fc_speed_qscore_2d_sink_opts cfg;
    char *output_root;
} speed_qscore_2d_sink;

// open the barcode/non-demux output file for the 2D table
static FILE *qts2d_open_output(const speed_qscore_2d_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (!dir) {
        FC_LOG_FATAL("Failed to create output directory for speed_qscore_2d (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "speed_qscore_2d.tsv");
    free(dir);
    ensure_parent_dir_exists(path);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open speed_qscore_2d output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

// write a materialized dense speed x qscore matrix to TSV
static void qts2d_write_matrix(const fc_speed_qscore_2d_matrix *matrix, FILE *fp)
{
    fprintf(fp, "speed_bin_start\tspeed_bin_end\tqscore_bin_start\tqscore_bin_end\tcount\n");

    double row_width = (matrix->x_upper - matrix->x_lower) / (double)matrix->x_bins;
    for (size_t speed_bin = 0; speed_bin < matrix->x_bins; speed_bin++) {
        double speed_start = matrix->x_lower + (double)speed_bin * row_width;
        double speed_end = (speed_bin + 1 == matrix->x_bins)
            ? matrix->x_upper
            : matrix->x_lower + (double)(speed_bin + 1) * row_width;
        for (size_t q = 0; q < matrix->qscore_bins; q++) {
            double q_start = (double)q * matrix->qscore_width;
            double q_end = q_start + matrix->qscore_width;
            double count = matrix->counts[q * matrix->x_bins + speed_bin];
            if (!isfinite(count) || count < 0.0) {
                count = 0.0;
            }
            fprintf(fp, "%.10g\t%.10g\t%.1f\t%.1f\t%.12g\n",
                speed_start, speed_end, q_start, q_end, count);
        }
    }
}

// materialize one barcode matrix and write dense rows to TSV
static void qts2d_write_barcode_output(const speed_qscore_2d_sink *sink, size_t barcode)
{
    fc_speed_qscore_2d_matrix matrix = { 0 };
    if (fc_speed_qscore_2d_materialize(sink->cfg.model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize speed_qscore_2d matrix for barcode %zu.", barcode);
    }

    FILE *fp = qts2d_open_output(sink, barcode);
    qts2d_write_matrix(&matrix, fp);
    fclose(fp);
    fc_speed_qscore_2d_matrix_destroy(&matrix);
}

// flush outputs for seen barcode states
static void speed_qscore_2d_sink_finish(void *ctx)
{
    speed_qscore_2d_sink *sink = (speed_qscore_2d_sink *)ctx;
    fc_speed_qscore_2d_skip_counts skipped = fc_speed_qscore_2d_skipped_counts(sink->cfg.model);
    if (skipped.duration_missing || skipped.duration_invalid) {
        FC_LOG_WARN(
            "Skipped %zu reads without du and %zu reads with invalid du for speed_qscore_2d.",
            skipped.duration_missing,
            skipped.duration_invalid);
    }

    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (fc_speed_qscore_2d_seen(sink->cfg.model, b)) {
                qts2d_write_barcode_output(sink, b);
            }
        }
    } else {
        qts2d_write_barcode_output(sink, 0);
    }
}

// destroy sink wrapper state without owning the shared model
static void speed_qscore_2d_sink_destroy(void *ctx)
{
    speed_qscore_2d_sink *sink = (speed_qscore_2d_sink *)ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_speed_qscore_2d_sink_create(
    const fc_speed_qscore_2d_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root || !opts->model) {
        return sink;
    }

    speed_qscore_2d_sink *state = fc_alloc(
        1,
        sizeof(*state),
        "speed_qscore_2d sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);

    sink.ctx = state;
    sink.finish = speed_qscore_2d_sink_finish;
    sink.destroy = speed_qscore_2d_sink_destroy;
    return sink;
}
