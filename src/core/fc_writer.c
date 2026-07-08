#include "core/fc_writer.h"

#include <stdio.h>
#include <stdlib.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_log.h"

fc_writer *fc_writer_create(const fc_writer_config *cfg)
{
    fc_writer *writer = fc_alloc(1, sizeof(*writer), "writer");
    const fc_stats_engine_config *stats_cfg = cfg ? cfg->stats : NULL;
    bool run_summary_enabled = cfg ? cfg->run_summary : false;
    bool channel_summary_enabled = cfg ? cfg->channel_summary : false;
    const fc_length_qscore_2d_config *length_qscore_cfg = cfg ? cfg->length_qscore : NULL;
    const fc_speed_qscore_2d_config *speed_qscore_cfg = cfg ? cfg->speed_qscore : NULL;
    const fc_time_bin_config *time_bin_cfg = cfg ? cfg->time_bin : NULL;
    bool use_barcodes = cfg ? cfg->use_barcodes : false;

    writer->use_barcodes = use_barcodes;
    writer->barcode_count = use_barcodes ? FC_MAX_BARCODES : 1;

    writer->stats_per_barcode = use_barcodes;
    writer->stats_cfg_set = (stats_cfg != NULL);
    if (writer->stats_cfg_set) {
        writer->stats_cfg = *stats_cfg;
        writer->stats_engine_count = writer->stats_per_barcode ? writer->barcode_count : 1;
        writer->stats_engines = fc_alloc(
            writer->stats_engine_count, sizeof(fc_stats_engine *), "stats engines");
        // Non-demux runs use a single shared engine; create it eagerly.
        if (!writer->stats_per_barcode && writer->stats_engine_count > 0) {
            writer->stats_engines[0] = fc_stats_engine_create(
                writer->stats_cfg.track_alignment,
                writer->stats_cfg.track_polya,
                writer->stats_cfg.polya_cover,
                writer->stats_cfg.polya_qual,
                writer->stats_cfg.polya_rev);
        }
    }
    if (run_summary_enabled) {
        writer->run_summary = fc_run_summary_create(writer->barcode_count);
        if (writer->run_summary == NULL) {
            FC_LOG_ERROR("Failed to create run-summary model.");
            fc_writer_destroy(writer);
            return NULL;
        }
    }
    if (length_qscore_cfg != NULL) {
        writer->length_qscore = fc_length_qscore_2d_create(
            writer->barcode_count,
            length_qscore_cfg->kll_k,
            length_qscore_cfg->qgrid_points,
            length_qscore_cfg->smoothing_points);
        if (writer->length_qscore == NULL) {
            FC_LOG_ERROR("Failed to create length/qscore 2D model.");
            fc_writer_destroy(writer);
            return NULL;
        }
    }
    if (speed_qscore_cfg != NULL) {
        writer->speed_qscore = fc_speed_qscore_2d_create(
            writer->barcode_count,
            speed_qscore_cfg->kll_k,
            speed_qscore_cfg->qgrid_points,
            speed_qscore_cfg->smoothing_points);
        if (writer->speed_qscore == NULL) {
            FC_LOG_ERROR("Failed to create qscore/translocation-speed 2D model.");
            fc_writer_destroy(writer);
            return NULL;
        }
    }
    if (channel_summary_enabled) {
        writer->channel_summary = fc_channel_summary_create(writer->barcode_count);
        if (writer->channel_summary == NULL) {
            FC_LOG_ERROR("Failed to create channel-summary model.");
            fc_writer_destroy(writer);
            return NULL;
        }
    }
    if (time_bin_cfg != NULL) {
        writer->time_bin = fc_time_bin_create(
            writer->barcode_count,
            time_bin_cfg->bin_size_s,
            time_bin_cfg->kll_k);
        if (writer->time_bin == NULL) {
            FC_LOG_ERROR("Failed to create time-bin model.");
            fc_writer_destroy(writer);
            return NULL;
        }
    }
    return writer;
}

void fc_writer_add_sink(fc_writer *writer, fc_sink sink)
{
    if (!writer) {
        return;
    }
    if (writer->sink_count == writer->sink_capacity) {
        size_t new_cap = writer->sink_capacity == 0 ? 4 : writer->sink_capacity * 2;
        writer->sinks = fc_realloc(writer->sinks, new_cap * sizeof(*writer->sinks), "writer sinks array");
        writer->sink_capacity = new_cap;
    }
    writer->sinks[writer->sink_count++] = sink;
}

void fc_writer_begin_file(fc_writer *writer, const char *fname)
{
    if (!writer) {
        return;
    }
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].begin_file) {
            writer->sinks[i].begin_file(writer->sinks[i].ctx, fname);
        }
    }
}

void fc_writer_add_read(fc_writer *writer, fc_read *read, const char *fname, const sam_hdr_t *hdr)
{
    if (!writer) {
        return;
    }
    fc_read_ctx ctx = { 0 };
    ctx.barcode = writer->use_barcodes ? read->meta->ibarcode : 0;
    if (writer->barcode_count > 0 && ctx.barcode >= writer->barcode_count) {
        FC_LOG_FATAL("read barcode (%zu) exceeds configured count (%zu)",
            ctx.barcode, writer->barcode_count);
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

    if (writer->length_qscore != NULL) {
        size_t barcode = writer->use_barcodes ? ctx.barcode : 0;
        fc_length_qscore_2d_add_read(
            writer->length_qscore,
            barcode,
            (double)read->bam->core.l_qseq,
            read->meta->mean_q);
    }

    if (writer->speed_qscore != NULL) {
        size_t barcode = writer->use_barcodes ? ctx.barcode : 0;
        fc_speed_qscore_2d_add_read(
            writer->speed_qscore,
            barcode,
            (double)read->bam->core.l_qseq,
            read->meta->duration_s,
            read->meta->mean_q);
    }

    if (writer->run_summary != NULL) {
        size_t barcode = writer->use_barcodes ? ctx.barcode : 0;
        fc_run_summary_add_read(writer->run_summary, barcode, read);
    }

    if (writer->channel_summary != NULL) {
        size_t barcode = writer->use_barcodes ? ctx.barcode : 0;
        fc_channel_summary_add_read(writer->channel_summary, barcode, read);
    }

    if (writer->time_bin != NULL) {
        size_t barcode = writer->use_barcodes ? ctx.barcode : 0;
        fc_time_bin_add_read(writer->time_bin, barcode, read);
    }

    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].add_read) {
            writer->sinks[i].add_read(writer->sinks[i].ctx, read, &ctx, fname, hdr);
        }
    }
}

void fc_writer_end_file(fc_writer *writer, const char *fname, const fc_file_stats *stats)
{
    if (!writer) {
        return;
    }
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].end_file) {
            writer->sinks[i].end_file(writer->sinks[i].ctx, fname, stats);
        }
    }
}

void fc_writer_finish(fc_writer *writer)
{
    if (!writer) {
        return;
    }
    for (size_t i = 0; i < writer->sink_count; i++) {
        if (writer->sinks[i].finish) {
            writer->sinks[i].finish(writer->sinks[i].ctx);
        }
    }
}

void fc_writer_destroy(fc_writer *writer)
{
    if (!writer) {
        return;
    }
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
    fc_run_summary_destroy(writer->run_summary);
    fc_channel_summary_destroy(writer->channel_summary);
    fc_length_qscore_2d_destroy(writer->length_qscore);
    fc_speed_qscore_2d_destroy(writer->speed_qscore);
    fc_time_bin_destroy(writer->time_bin);
    free(writer->sinks);
    free(writer);
}
