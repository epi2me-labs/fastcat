#include "fc_sinks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "sink/common.h"
#include "htslib/sam.h"
#include "../fc_log.h"

typedef struct {
    fc_perread_sink_opts cfg;  // cached configuration
    FILE **fp;                 // one per barcode slot
    bool *header_written;
    char *sample;              // duplicated if provided
    char *output_root;         // top-level output root (stripped)
} perread_sink;

static void perread_sink_write_header(FILE *fp, const char *sample, bool include_alignment) {
    if (!fp) return;
    fprintf(fp, "read_id\tfilename\trunid\t");
    if (sample) fprintf(fp, "sample_name\t");
    fprintf(fp, "read_length\tmean_quality\tchannel\tread_number\tstart_time");
    if (include_alignment) {
        fprintf(fp,
            "\tref\tcoverage\tref_coverage\tqstart\tqend\t"
            "rstart\trend\taligned_ref_len\tdirection\tlength\t"
            "match\tins\tdel\tsub\tiden\tacc\tduplex");
    }
    fprintf(fp, "\n");
}

static FILE *perread_sink_ensure_fp(perread_sink *sink, size_t barcode) {
    if (barcode >= sink->cfg.barcode_count) return NULL;
    if (sink->fp[barcode]) return sink->fp[barcode];

    const char *leaf = "per_read.tsv";
    char *dir = sink->cfg.use_barcodes
        ? fc_sink_ensure_barcode_dir(sink->output_root, barcode)
        : fc_strip_path(sink->output_root);
    if (!dir) {
        FC_LOG_ERROR("Failed to create barcode directory for per-read output (barcode %zu).", barcode);
        exit(EXIT_FAILURE);
    }
    char *path = fc_join_path(dir, leaf);
    free(dir);
    ensure_parent_dir_exists(path);
    sink->fp[barcode] = fopen(path, "w");
    if (!sink->fp[barcode]) {
        FC_LOG_ERROR("Cannot open per-read file '%s' for writing.", path);
        free(path);
        exit(EXIT_FAILURE);
    }
    free(path);
    perread_sink_write_header(sink->fp[barcode], sink->sample, sink->cfg.include_alignment);
    sink->header_written[barcode] = true;
    return sink->fp[barcode];
}

static void perread_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr) {
    perread_sink *sink = (perread_sink *)ctx;
    size_t barcode = sink->cfg.use_barcodes && read_ctx ? read_ctx->barcode : 0;
    FILE *out = perread_sink_ensure_fp(sink, barcode);
    if (!out) return;
    if (sink->cfg.include_alignment) {
        if (!read_ctx) return;
        if (read_ctx->stats_kind != FC_STATS_ALIGNED && read_ctx->stats_kind != FC_STATS_UNMAPPED) {
            return;
        }
    }

    const char *sample = sink->sample ? sink->sample : "";
    const char *sample_tab = sink->sample ? "\t" : "";
    const char *runid = read->meta->runid ? read->meta->runid : "";
    const char *start_time = read->meta->start_time ? read->meta->start_time : "";
    uint32_t read_length = read->bam->core.l_qseq;
    float mean_quality = read->meta->mean_q;

    fprintf(out, "%s\t%s\t%s%s%s\t%u\t%.2f\t%zu\t%zu\t%s",
        bam_get_qname(read->bam), fname, runid,
        sample_tab, sample, read_length,
        mean_quality, read->meta->channel, read->meta->read_number, start_time);

    if (!sink->cfg.include_alignment) {
        fprintf(out, "\n");
        return;
    }

    const char *ref = "*";
    if (hdr && read->bam->core.tid >= 0) {
        ref = sam_hdr_tid2name(hdr, read->bam->core.tid);
    }

    if (read_ctx->stats_kind == FC_STATS_UNMAPPED) {
        fprintf(out,
            "\t%s\tnan\tnan\tnan\tnan\tnan\tnan\t0\t*\t0\t0\t0\t0\t0\tnan\tnan\t%d\n",
            ref, read->meta->duplex_code);
        return;
    }

    fprintf(out,
        "\t%s\t%.4f\t%.4f\t%zu\t%zu\t%zu\t%zu\t%zu\t%c\t%zu\t%zu\t%zu\t%zu\t%zu\t%.2f\t%.2f\t%d\n",
        ref,
        read_ctx->aligned.coverage, read_ctx->aligned.ref_cover,
        read_ctx->aligned.qstart, read_ctx->aligned.qend,
        read_ctx->aligned.rstart, read_ctx->aligned.rend,
        read_ctx->aligned.aligned_ref_len, read_ctx->aligned.direction,
        read_ctx->aligned.length,
        read_ctx->aligned.match, read_ctx->aligned.ins,
        read_ctx->aligned.delt, read_ctx->aligned.sub,
        read_ctx->aligned.iden, read_ctx->aligned.acc,
        read->meta->duplex_code);
}

static void perread_sink_destroy(void *ctx) {
    perread_sink *sink = (perread_sink *)ctx;
    if (sink->fp) {
        for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
            if (sink->fp[i]) fclose(sink->fp[i]);
        }
        free(sink->fp);
    }
    free(sink->header_written);
    free(sink->sample);
    free(sink->output_root);
    free(sink);
}

fc_sink fc_perread_sink_create(const fc_perread_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts || !opts->output_root) return sink;

    perread_sink *state = calloc(1, sizeof(*state));
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output_root = fc_strip_path(opts->output_root);
    state->fp = calloc(state->cfg.barcode_count, sizeof(FILE *));
    state->header_written = calloc(state->cfg.barcode_count, sizeof(bool));
    if (opts->sample) {
        state->sample = strdup(opts->sample);
    }

    sink.ctx = state;
    sink.add_read = perread_sink_add_read;
    sink.destroy = perread_sink_destroy;
    return sink;
}
