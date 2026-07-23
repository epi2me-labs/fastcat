#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"
#include "core/fc_log.h"
#include "core/fc_sinks.h"
#include "sink/common.h"

typedef struct flag_stats {
    size_t n_refs;
    size_t **counts;
    size_t *unmapped;
} flag_stats;

// Allocate and destroy flag statistics structures.
static flag_stats *create_flag_stats(size_t n_refs, bool store_unmapped)
{
    flag_stats *stats = fc_alloc(1, sizeof(flag_stats), "flagstat");
    stats->n_refs = n_refs;
    // For unaligned-only BAMs (n_refs=0), we only track unmapped stats
    stats->counts = n_refs > 0 ? fc_alloc(n_refs, sizeof(size_t *), "flagstat") : NULL;
    stats->unmapped = store_unmapped ? fc_alloc(FC_FLAGSTAT_COUNT, sizeof(size_t), "flagstat") : NULL;

    for (size_t i = 0; i < n_refs; i++) {
        stats->counts[i] = fc_alloc(FC_FLAGSTAT_COUNT, sizeof(size_t), "flagstat");
    }
    return stats;
}

static void destroy_flag_stats(flag_stats *stats)
{
    if (!stats) {
        return;
    }
    for (size_t i = 0; i < stats->n_refs; i++) {
        free(stats->counts[i]);
    }
    free(stats->counts);
    free(stats->unmapped);
    free(stats);
}

static const size_t flag_masks[FC_FLAGSTAT_COUNT] = {
    [FC_FLAGSTAT_SECONDARY] = BAM_FSECONDARY,
    [FC_FLAGSTAT_SUPPLEMENTARY] = BAM_FSUPPLEMENTARY,
    [FC_FLAGSTAT_UNMAPPED] = BAM_FUNMAP,
    [FC_FLAGSTAT_QCFAIL] = BAM_FQCFAIL,
    [FC_FLAGSTAT_DUPLICATE] = BAM_FDUP
};

static inline void flagstats_add_counts(const bam1_t *b, size_t *counts, const int duplex_code)
{
    counts[FC_FLAGSTAT_TOTAL] += 1;
    counts[FC_FLAGSTAT_PRIMARY] += fc_is_primary(b);
    for (size_t i = FC_FLAGSTAT_SECONDARY; i <= FC_FLAGSTAT_DUPLICATE; ++i) {
        counts[i] += ((b->core.flag & flag_masks[i]) != 0);
    }
    counts[FC_FLAGSTAT_DUPLEX] += (duplex_code == 1);
    counts[FC_FLAGSTAT_DUPLEX_FORMING] += (duplex_code == -1);
}

typedef struct {
    fc_flagstats_sink_opts cfg; // cached configuration
    char *sample_dup;
    char *region_label_dup;
    flag_stats *global_stats;
    flag_stats **bc_stats;
    FILE *global_fp;
    FILE **bc_fp;
    char *output_root;
    size_t n_refs;
} flagstats_sink;

static flag_stats *ensure_flag_stats(flagstats_sink *sink, size_t barcode)
{
    if (barcode >= sink->cfg.barcode_count) {
        return NULL;
    }
    if (!sink->bc_stats[barcode]) {
        sink->bc_stats[barcode] = create_flag_stats(sink->n_refs, true);
    }
    return sink->bc_stats[barcode];
}

static void flagstats_sink_add_read(void *ctx, const fc_read *fr, const fc_read_ctx *read_ctx,
    const char *fname, const sam_hdr_t *hdr)
{
    flagstats_sink *sink = (flagstats_sink *)ctx;
    (void)fname;
    (void)hdr;
    size_t barcode = (sink->cfg.use_barcodes && read_ctx) ? read_ctx->barcode : 0;

    const bam1_t *b = fr->bam;
    bool is_unmapped = fc_is_unmapped(b);
    // update global
    if (!is_unmapped) {
        size_t tid = (size_t)b->core.tid;
        if (!sink->cfg.region_mode && (tid >= sink->n_refs)) {
            return;
        }
        size_t *counts = sink->cfg.region_mode ? sink->global_stats->counts[0] : sink->global_stats->counts[tid];
        flagstats_add_counts(b, counts, fr->meta->duplex_code);
    }
    if (sink->global_stats->unmapped && is_unmapped) {
        flagstats_add_counts(b, sink->global_stats->unmapped, fr->meta->duplex_code);
    }

    if (sink->cfg.use_barcodes && barcode < sink->cfg.barcode_count) {
        flag_stats *bstats = ensure_flag_stats(sink, barcode);
        if (bstats) {
            if (!is_unmapped) {
                size_t tid = (size_t)b->core.tid;
                if (!sink->cfg.region_mode && (tid < sink->n_refs)) {
                    size_t *counts = sink->cfg.region_mode ? bstats->counts[0] : bstats->counts[tid];
                    flagstats_add_counts(b, counts, fr->meta->duplex_code);
                }
            }
            if (bstats->unmapped && is_unmapped) {
                flagstats_add_counts(b, bstats->unmapped, fr->meta->duplex_code);
            }
        }
    }
}

static void flagstats_write_header(FILE *fh, const char *sample)
{
    if (!fh) {
        return;
    }
    fprintf(fh, "ref%s", sample ? "\tsample_name" : "");
    for (size_t i = 0; i < FC_FLAGSTAT_COUNT; ++i) {
        fprintf(fh, "\t%s", fc_flagstat_column_names[i]);
    }
    fputc('\n', fh);
}

static void flagstats_write_line(FILE *fh, const char *sample, const char *label, size_t *stats)
{
    if (!fh) {
        return;
    }
    fprintf(fh, "%s", label);
    if (sample) {
        fprintf(fh, "\t%s", sample);
    }
    for (size_t i = 0; i < FC_FLAGSTAT_COUNT; ++i) {
        fprintf(fh, "\t%zu", stats[i]);
    }
    fputc('\n', fh);
}

static void flagstats_write_all(flagstats_sink *sink, flag_stats *stats, FILE *fp)
{
    if (!stats || !fp) {
        return;
    }

    if (sink->cfg.region_mode) {
        const char *label = sink->cfg.region_label ? sink->cfg.region_label : "region";
        flagstats_write_line(fp, sink->cfg.sample, label, stats->counts[0]);
    } else {
        for (size_t i = 0; i < stats->n_refs; ++i) {
            const char *label = NULL;
            if (sink->cfg.hdr && i < (size_t)sink->cfg.hdr->n_targets) {
                label = sam_hdr_tid2name((sam_hdr_t *)sink->cfg.hdr, (int)i);
            }
            char tmp[32];
            if (!label) {
                FC_SNPRINTF_CHECK(tmp, sizeof(tmp), "%zu", i);
                label = tmp;
            }
            flagstats_write_line(fp, sink->cfg.sample, label, stats->counts[i]);
        }
    }
    if (stats->unmapped) {
        flagstats_write_line(fp, sink->cfg.sample, "*", stats->unmapped);
    }
}

static FILE *flagstats_open_file(flagstats_sink *sink, const char *path)
{
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_FATAL("Cannot open file '%s' for writing.", path);
    }
    flagstats_write_header(fp, sink->cfg.sample);
    return fp;
}

static void flagstats_sink_finish(void *ctx)
{
    flagstats_sink *s = (flagstats_sink *)ctx;
    if (!s) {
        return;
    }
    if (s->global_fp && s->global_stats) {
        flagstats_write_all(s, s->global_stats, s->global_fp);
    }
    if (s->cfg.use_barcodes) {
        const char *leaf = "flagstats.tsv";
        for (size_t b = 0; b < s->cfg.barcode_count; ++b) {
            if (!s->bc_stats[b]) {
                continue;
            }
            if (!s->bc_fp[b]) {
                char *dir = fc_sink_ensure_barcode_dir(s->output_root, b);
                if (!dir) {
                    FC_LOG_FATAL("Failed to create barcode directory for flagstats (barcode %zu).", b);
                }
                char *path = fc_join_path(dir, leaf);
                free(dir);
                s->bc_fp[b] = flagstats_open_file(s, path);
                free(path);
            }
            flagstats_write_all(s, s->bc_stats[b], s->bc_fp[b]);
        }
    }
}

static void flagstats_sink_destroy(void *ctx)
{
    flagstats_sink *sink = (flagstats_sink *)ctx;
    if (sink->global_fp) {
        fclose(sink->global_fp);
    }
    if (sink->bc_fp) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->bc_fp[i]) {
                fclose(sink->bc_fp[i]);
            }
        }
    }
    free(sink->sample_dup);
    free(sink->region_label_dup);
    destroy_flag_stats(sink->global_stats);
    if (sink->bc_stats) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            destroy_flag_stats(sink->bc_stats[i]);
        }
    }
    free(sink->bc_stats);
    free(sink->bc_fp);
    free(sink->output_root);
    free(sink);
}

fc_sink fc_flagstats_sink_create(const fc_flagstats_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts || !opts->output_root) {
        return sink;
    }

    flagstats_sink *state = fc_alloc(1, sizeof(*state), "flagstats sink");
    state->cfg = *opts;
    state->cfg.sample = NULL;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->n_refs = state->cfg.region_mode ? 1 : (state->cfg.hdr ? (size_t)state->cfg.hdr->n_targets : 0);
    state->global_stats = create_flag_stats(state->n_refs, true);
    state->bc_stats = fc_alloc(state->cfg.barcode_count, sizeof(flag_stats *), "flagstats barcode stats");
    state->bc_fp = fc_alloc(state->cfg.barcode_count, sizeof(FILE *), "flagstats barcode files");
    state->output_root = fc_strip_path(opts->output_root);
    if (opts->sample) {
        state->sample_dup = fc_strdup(opts->sample, "sample name for flagstats sink");
        state->cfg.sample = state->sample_dup;
    }
    if (opts->region_label) {
        state->region_label_dup = fc_strdup(opts->region_label, "region label for flagstats sink");
        state->cfg.region_label = state->region_label_dup;
    }
    if (state->output_root) {
        char *path = fc_join_path(state->output_root, "flagstats.tsv");
        state->global_fp = flagstats_open_file(state, path);
        free(path);
    }

    sink.ctx = state;
    sink.add_read = flagstats_sink_add_read;
    sink.finish = flagstats_sink_finish;
    sink.destroy = flagstats_sink_destroy;
    return sink;
}
