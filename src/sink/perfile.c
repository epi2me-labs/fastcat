#include "fc_sinks.h"

#include <ctype.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "fastcat/parsing.h"
#include "sink/common.h"

typedef struct {
    FILE *fp;
    char *sample;
} perfile_sink;

static void perfile_sink_end_file(void *ctx, const char *fname, const fc_file_stats *stats) {
    perfile_sink *sink = (perfile_sink *)ctx;
    if (!sink->fp || !stats) return;

    fprintf(sink->fp, "%s\t", fname);
    if (sink->sample) fprintf(sink->fp, "%s\t", sink->sample);
    if (stats->n_seqs == 0) {
        fprintf(sink->fp, "0\t0\t0\t0\t0.00");
    } else {
        fprintf(sink->fp, "%zu\t%zu\t%zu\t%zu\t%.2f",
            stats->n_seqs, stats->n_bases, stats->min_length, stats->max_length, stats->mean_quality);
    }
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        fprintf(sink->fp, "\t%" PRIu64, stats->failures[i]);
    }
    fprintf(sink->fp, "\n");
}

static void perfile_sink_destroy(void *ctx) {
    perfile_sink *sink = (perfile_sink *)ctx;
    if (sink->fp) fclose(sink->fp);
    free(sink->sample);
    free(sink);
}

fc_sink fc_perfile_sink_create(const char *output_root, const char *sample) {
    fc_sink sink = {0};
    if (!output_root) return sink;

    perfile_sink *state = calloc(1, sizeof(*state));
    char *path = fc_sink_join_path(output_root, "per_file.tsv");
    ensure_parent_dir_exists(path);
    state->fp = fopen(path, "w");
    if (state->fp == NULL) {
        fprintf(stderr, "Error opening per-file file '%s' for writing.\n", path);
        exit(EXIT_FAILURE);
    }
    free(path);
    if (sample && sample[0]) {
        state->sample = strdup(sample);
    }

    fprintf(state->fp, "filename\t");
    if (state->sample) fprintf(state->fp, "sample_name\t");
    fprintf(state->fp, "n_seqs\tn_bases\tmin_length\tmax_length\tmean_quality");
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        const char *src = failure_type[i];
        size_t len = strlen(src);
        char *buffer = calloc(len + 1, sizeof(char));
        for (size_t j = 0; j < len; ++j) {
            buffer[j] = tolower((unsigned char)src[j]);
        }
        fprintf(state->fp, "\t%s", buffer);
        free(buffer);
    }
    fprintf(state->fp, "\n");

    sink.ctx = state;
    sink.end_file = perfile_sink_end_file;
    sink.destroy = perfile_sink_destroy;
    return sink;
}
