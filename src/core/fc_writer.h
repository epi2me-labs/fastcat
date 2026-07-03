#ifndef FASTCAT_FC_WRITER_H
#define FASTCAT_FC_WRITER_H

#include <stddef.h>

#include "core/fc_iter.h"
#include "core/fc_length_qscore_2d.h"
#include "core/fc_stats_engine.h"
#include "core/fc_time_bin.h"

// NOTE: The fc_writer, fc_sink, and fc_stats_engine interfaces are designed
// for single-threaded record processing only. They are NOT thread-safe.
// Only the underlying BAM writer (via htslib thread pool) may use multiple threads.

typedef struct {
    size_t barcode;
    bool is_unmapped;
    bool is_primary_qc;
    fc_stats_read_kind stats_kind;
    fc_stats_aligned aligned;
} fc_read_ctx;

typedef struct {
    void *ctx;
    void (*begin_file)(void *ctx, const char *fname);
    // `read` is only valid for the duration of this callback. Sinks must not retain
    // pointers into `read` (or `read->meta`) after returning.
    void (*add_read)(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr);
    void (*end_file)(void *ctx, const char *fname, const fc_file_stats *stats);
    void (*finish)(void *ctx);
    void (*destroy)(void *ctx);
} fc_sink;

typedef struct {
    bool track_alignment;
    bool track_polya;
    float polya_cover;
    float polya_qual;
    bool polya_rev;
} fc_stats_engine_config;

typedef struct {
    bool use_barcodes;
    const fc_stats_engine_config *stats;
    const fc_length_qscore_2d_config *length_qscore;
    const fc_time_bin_config *time_bin;
} fc_writer_config;

typedef struct {
    fc_sink *sinks;
    size_t sink_count;
    size_t sink_capacity;
    bool use_barcodes;
    size_t barcode_count;
    bool stats_per_barcode;
    bool stats_cfg_set;
    fc_stats_engine **stats_engines;
    size_t stats_engine_count;
    fc_stats_engine_config stats_cfg;
    fc_length_qscore_2d *length_qscore;
    fc_time_bin *time_bin;
} fc_writer;

fc_writer *fc_writer_create(const fc_writer_config *cfg);
void fc_writer_add_sink(fc_writer *writer, fc_sink sink);

void fc_writer_begin_file(fc_writer *writer, const char *fname);
// Dispatches a single record to all sinks. `read` must remain valid for this call only.
void fc_writer_add_read(fc_writer *writer, fc_read *read, const char *fname, const sam_hdr_t *hdr);
void fc_writer_end_file(fc_writer *writer, const char *fname, const fc_file_stats *stats);
void fc_writer_finish(fc_writer *writer);
void fc_writer_destroy(fc_writer *writer);

#endif
