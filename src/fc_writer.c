#include "fc_writer.h"

#include <stdio.h>
#include <stdlib.h>

#include "fc_filter.h"

fc_writer *fc_writer_create(void) {
    fc_writer *writer = calloc(1, sizeof(*writer));
    writer->barcode_count = 1;
    return writer;
}

void fc_writer_add_sink(fc_writer *writer, fc_sink sink) {
    if (!writer) return;
    if (writer->sink_count == writer->sink_capacity) {
        size_t new_cap = writer->sink_capacity == 0 ? 4 : writer->sink_capacity * 2;
        fc_sink *new_sinks = realloc(writer->sinks, new_cap * sizeof(*new_sinks));
        if (!new_sinks) {
            exit(EXIT_FAILURE);
        }
        writer->sinks = new_sinks;
        writer->sink_capacity = new_cap;
    }
    writer->sinks[writer->sink_count++] = sink;
}

void fc_writer_set_barcode_count(fc_writer *writer, size_t barcode_count, bool use_barcodes) {
    if (!writer) return;
    writer->barcode_count = barcode_count;
    writer->use_barcodes = use_barcodes;
}

void fc_writer_configure_stats(fc_writer *writer, const fc_stats_engine_config *config, bool per_barcode) {
    if (!writer || !config) return;
    if (writer->stats_engines) {
        for (size_t i = 0; i < writer->stats_engine_count; ++i) {
            if (writer->stats_engines[i] != NULL) {
                fc_stats_engine_destroy(writer->stats_engines[i]);
            }
        }
        free(writer->stats_engines);
    }
    writer->stats_per_barcode = per_barcode;
    writer->stats_cfg_set = true;
    writer->stats_cfg = *config;
    writer->stats_engine_count = per_barcode ? writer->barcode_count : 1;
    writer->stats_engines = calloc(writer->stats_engine_count, sizeof(*writer->stats_engines));
    if (!writer->stats_engines) {
        fprintf(stderr, "ERROR: failed to allocate stats engines\n");
        exit(EXIT_FAILURE);
    }
    if (config->create_default && writer->stats_engine_count > 0) {
        writer->stats_engines[0] = fc_stats_engine_create(
            config->track_unmapped,
            config->track_alignment,
            config->track_polya,
            config->polya_cover,
            config->polya_qual,
            config->polya_rev);
    }
}

fc_stats_engine **fc_writer_get_stats_engines(fc_writer *writer, size_t *engine_count) {
    if (!writer) return NULL;
    if (engine_count) *engine_count = writer->stats_engine_count;
    return writer->stats_engines;
}

void fc_writer_begin_file(fc_writer *writer, const char *fname) {
    if (!writer) return;
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].begin_file) {
            writer->sinks[i].begin_file(writer->sinks[i].ctx, fname);
        }
    }
}

void fc_writer_add_read(fc_writer *writer, fc_read *read, const char *fname, const sam_hdr_t *hdr) {
    if (!writer) return;
    fc_read_ctx ctx = {0};
    ctx.barcode = writer->use_barcodes ? read->meta->ibarcode : 0;
    if (writer->barcode_count > 0 && ctx.barcode >= writer->barcode_count) {
        fprintf(stderr, "ERROR: read barcode (%zu) exceeds configured count (%zu)\n",
            ctx.barcode, writer->barcode_count);
        exit(EXIT_FAILURE);
    }
    ctx.is_unmapped = fc_is_unmapped(read->bam);
    ctx.is_primary_qc = fc_is_primary_qc_pass(read->bam);
    ctx.stats_kind = FC_STATS_SKIP;

    if (writer->stats_engines) {
        size_t engine_idx = writer->stats_per_barcode ? ctx.barcode : 0;
        if (engine_idx < writer->stats_engine_count) {
            fc_stats_engine *engine = writer->stats_engines[engine_idx];
            if (engine == NULL && writer->stats_cfg_set) {
                engine = fc_stats_engine_create(
                    writer->stats_cfg.track_unmapped,
                    writer->stats_cfg.track_alignment,
                    writer->stats_cfg.track_polya,
                    writer->stats_cfg.polya_cover,
                    writer->stats_cfg.polya_qual,
                    writer->stats_cfg.polya_rev);
                writer->stats_engines[engine_idx] = engine;
            }
            if (engine != NULL) {
                ctx.stats_kind = fc_stats_engine_add_read(engine, read, hdr, &ctx.aligned);
            }
        }
    }

    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].add_read) {
            writer->sinks[i].add_read(writer->sinks[i].ctx, read, &ctx, fname, hdr);
        }
    }
}

void fc_writer_end_file(fc_writer *writer, const char *fname, const fc_file_stats *stats) {
    if (!writer) return;
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].end_file) {
            writer->sinks[i].end_file(writer->sinks[i].ctx, fname, stats);
        }
    }
}

void fc_writer_finish(fc_writer *writer) {
    if (!writer) return;
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].finish) {
            writer->sinks[i].finish(writer->sinks[i].ctx);
        }
    }
}

void fc_writer_destroy(fc_writer *writer) {
    if (!writer) return;
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].destroy) {
            writer->sinks[i].destroy(writer->sinks[i].ctx);
        }
    }
    if (writer->stats_engines) {
        for (size_t i = 0; i < writer->stats_engine_count; ++i) {
            if (writer->stats_engines[i] != NULL) {
                fc_stats_engine_destroy(writer->stats_engines[i]);
            }
        }
        free(writer->stats_engines);
    }
    free(writer->sinks);
    free(writer);
}
