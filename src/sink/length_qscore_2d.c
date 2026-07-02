#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/fc_common.h"
#include "core/fc_length_qscore_2d.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct {
    fc_length_qscore_2d_sink_opts cfg;
    char *output_root;
} length_qscore_2d_sink;

// open the barcode/non-demux output file for the 2D table
static FILE *lq2d_open_output(const length_qscore_2d_sink *sink, size_t barcode)
{
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (!dir) {
        FC_LOG_FATAL("Failed to create output directory for length_qscore_2d (barcode %zu).", barcode);
    }

    char *path = fc_join_path(dir, "length_qscore_2d.tsv");
    free(dir);
    ensure_parent_dir_exists(path);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open length_qscore_2d output file '%s' for writing.", path);
    }
    free(path);
    return fp;
}

// write a materialized dense length x qscore matrix to TSV
static void lq2d_write_matrix(const fc_length_qscore_2d_matrix *matrix, FILE *fp)
{
    fprintf(fp, "length_bin_start\tlength_bin_end\tqscore_bin_start\tqscore_bin_end\tcount\n");

    double row_width = matrix->length_upper / (double)matrix->length_bins;
    for (size_t len_bin = 0; len_bin < matrix->length_bins; len_bin++) {
        double len_start = (double)len_bin * row_width;
        double len_end = (len_bin + 1 == matrix->length_bins)
            ? matrix->length_upper
            : (double)(len_bin + 1) * row_width;
        for (size_t q = 0; q < matrix->qscore_bins; q++) {
            double q_start = (double)q * matrix->qscore_width;
            double q_end = q_start + matrix->qscore_width;
            double count = matrix->counts[q * matrix->length_bins + len_bin];
            if (!isfinite(count) || count < 0.0) {
                count = 0.0;
            }
            fprintf(fp, "%.10g\t%.10g\t%.1f\t%.1f\t%.12g\n", len_start, len_end, q_start, q_end, count);
        }
    }
}

// materialize one barcode matrix and write dense rows to TSV
static void lq2d_write_barcode_output(const length_qscore_2d_sink *sink, size_t barcode)
{
    fc_length_qscore_2d_matrix matrix = { 0 };
    if (fc_length_qscore_2d_materialize(sink->cfg.model, barcode, &matrix) != 0) {
        FC_LOG_FATAL("Failed to materialize length_qscore_2d matrix for barcode %zu.", barcode);
    }

    FILE *fp = lq2d_open_output(sink, barcode);
    lq2d_write_matrix(&matrix, fp);
    fclose(fp);
    fc_length_qscore_2d_matrix_destroy(&matrix);
}

// flush outputs for active barcode states
static void length_qscore_2d_sink_finish(void *ctx)
{
    length_qscore_2d_sink *sink = (length_qscore_2d_sink *)ctx;
    if (sink->cfg.use_barcodes) {
        for (size_t b = 0; b < sink->cfg.barcode_count; b++) {
            if (fc_length_qscore_2d_active(sink->cfg.model, b)) {
                lq2d_write_barcode_output(sink, b);
            }
        }
    } else {
        lq2d_write_barcode_output(sink, 0);
    }
}

// destroy sink wrapper state without owning the shared model
static void length_qscore_2d_sink_destroy(void *ctx)
{
    length_qscore_2d_sink *sink = (length_qscore_2d_sink *)ctx;
    free(sink->output_root);
    free(sink);
}

fc_sink fc_length_qscore_2d_sink_create(const fc_length_qscore_2d_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root || !opts->model) {
        return sink;
    }

    length_qscore_2d_sink *state = fc_alloc(1, sizeof(*state), "length_qscore_2d sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);

    sink.ctx = state;
    sink.finish = length_qscore_2d_sink_finish;
    sink.destroy = length_qscore_2d_sink_destroy;
    return sink;
}
