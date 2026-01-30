#include "fc_sinks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "kh_counter.h"
#include "sink/common.h"
#include "../fc_log.h"

typedef struct {
    FILE *runids;
    FILE *basecallers;
    char *sample;
    kh_counter_t *run_ids;
    kh_counter_t *basecallers_counter;
} counter_sink;

static void counter_reset(counter_sink *sink) {
    if (sink->run_ids) kh_counter_destroy(sink->run_ids);
    if (sink->basecallers_counter) kh_counter_destroy(sink->basecallers_counter);
    sink->run_ids = kh_counter_init();
    sink->basecallers_counter = kh_counter_init();
}

static void counter_sink_begin_file(void *ctx, const char *fname) {
    counter_sink *sink = (counter_sink *)ctx;
    (void)fname;
    counter_reset(sink);
}

static void counter_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr) {
    counter_sink *sink = (counter_sink *)ctx;
    (void)fname;
    (void)hdr;
    (void)read_ctx;
    if (!sink->runids && !sink->basecallers) return;

    const char *runid = read->meta->runid ? read->meta->runid : "";
    const char *basecaller = read->meta->basecaller ? read->meta->basecaller : "";
    if (sink->runids) {
        kh_counter_increment(sink->run_ids, (char *)runid);
    }
    if (sink->basecallers) {
        kh_counter_increment(sink->basecallers_counter, (char *)basecaller);
    }
}

static void counter_sink_end_file(void *ctx, const char *fname, const fc_file_stats *stats) {
    counter_sink *sink = (counter_sink *)ctx;
    (void)stats;
    if (sink->runids) {
        for (khiter_t k = 0; k < kh_end(sink->run_ids); ++k) {
            if (kh_exist(sink->run_ids, k)) {
                fprintf(sink->runids, "%s\t", fname);
                if (sink->sample) fprintf(sink->runids, "%s\t", sink->sample);
                fprintf(sink->runids, "%s\t%d\n", kh_key(sink->run_ids, k), kh_val(sink->run_ids, k));
            }
        }
    }
    if (sink->basecallers) {
        for (khiter_t k = 0; k < kh_end(sink->basecallers_counter); ++k) {
            if (kh_exist(sink->basecallers_counter, k)) {
                fprintf(sink->basecallers, "%s\t", fname);
                if (sink->sample) fprintf(sink->basecallers, "%s\t", sink->sample);
                fprintf(sink->basecallers, "%s\t%d\n", kh_key(sink->basecallers_counter, k), kh_val(sink->basecallers_counter, k));
            }
        }
    }
}

static void counter_sink_destroy(void *ctx) {
    counter_sink *sink = (counter_sink *)ctx;
    if (sink->runids) fclose(sink->runids);
    if (sink->basecallers) fclose(sink->basecallers);
    if (sink->run_ids) kh_counter_destroy(sink->run_ids);
    if (sink->basecallers_counter) kh_counter_destroy(sink->basecallers_counter);
    free(sink->sample);
    free(sink);
}

fc_sink fc_counter_sink_create(const char *output_root, const char *sample) {
    fc_sink sink = {0};
    if (!output_root) return sink;

    counter_sink *state = calloc(1, sizeof(*state));
    if (sample && sample[0]) {
        state->sample = strdup(sample);
    }
    char *runids_path = fc_sink_join_path(output_root, "runids.tsv");
    ensure_parent_dir_exists(runids_path);
    state->runids = fopen(runids_path, "w");
    if (state->runids == NULL) {
        FC_LOG_ERROR("Error opening runids file '%s' for writing.", runids_path);
        exit(EXIT_FAILURE);
    }
    fprintf(state->runids, "filename\t");
    if (state->sample) fprintf(state->runids, "sample_name\t");
    fprintf(state->runids, "run_id\tcount\n");
    free(runids_path);

    char *basecallers_path = fc_sink_join_path(output_root, "basecallers.tsv");
    ensure_parent_dir_exists(basecallers_path);
    state->basecallers = fopen(basecallers_path, "w");
    if (state->basecallers == NULL) {
        FC_LOG_ERROR("Error opening basecallers file '%s' for writing.", basecallers_path);
        exit(EXIT_FAILURE);
    }
    fprintf(state->basecallers, "filename\t");
    if (state->sample) fprintf(state->basecallers, "sample_name\t");
    fprintf(state->basecallers, "basecaller\tcount\n");
    free(basecallers_path);

    state->run_ids = kh_counter_init();
    state->basecallers_counter = kh_counter_init();

    sink.ctx = state;
    sink.begin_file = counter_sink_begin_file;
    sink.add_read = counter_sink_add_read;
    sink.end_file = counter_sink_end_file;
    sink.destroy = counter_sink_destroy;
    return sink;
}
