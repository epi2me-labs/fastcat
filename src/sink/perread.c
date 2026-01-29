#include "fc_sinks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "sink/common.h"
#include "htslib/sam.h"

typedef struct {
    FILE *fp;
    char *sample;
    bool include_alignment;
    bool output_unmapped;
} perread_sink;

static void perread_sink_write_header(perread_sink *sink) {
    if (!sink || !sink->fp) return;
    fprintf(sink->fp, "read_id\tfilename\trunid\t");
    if (sink->sample) fprintf(sink->fp, "sample_name\t");
    fprintf(sink->fp, "read_length\tmean_quality\tchannel\tread_number\tstart_time");
    if (sink->include_alignment) {
        fprintf(sink->fp,
            "\tref\tcoverage\tref_coverage\tqstart\tqend\t"
            "rstart\trend\taligned_ref_len\tdirection\tlength\t"
            "match\tins\tdel\tsub\tiden\tacc\tduplex");
    }
    fprintf(sink->fp, "\n");
}

static void perread_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr) {
    perread_sink *sink = (perread_sink *)ctx;
    if (!sink->fp) return;
    if (sink->include_alignment) {
        if (!read_ctx) return;
        if (read_ctx->stats_kind == FC_STATS_UNMAPPED && !sink->output_unmapped) {
            return;
        }
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

    fprintf(sink->fp, "%s\t%s\t%s%s%s\t%u\t%.2f\t%zu\t%zu\t%s",
        bam_get_qname(read->bam), fname, runid,
        sample_tab, sample, read_length,
        mean_quality, read->meta->channel, read->meta->read_number, start_time);

    if (!sink->include_alignment) {
        fprintf(sink->fp, "\n");
        return;
    }

    const char *ref = "*";
    if (hdr && read->bam->core.tid >= 0) {
        ref = sam_hdr_tid2name(hdr, read->bam->core.tid);
    }

    if (read_ctx->stats_kind == FC_STATS_UNMAPPED) {
        fprintf(sink->fp,
            "\t%s\tnan\tnan\tnan\tnan\tnan\tnan\t0\t*\t0\t0\t0\t0\t0\tnan\tnan\t%d\n",
            ref, read->meta->duplex_code);
        return;
    }

    fprintf(sink->fp,
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
    if (sink->fp) fclose(sink->fp);
    free(sink->sample);
    free(sink);
}

fc_sink fc_perread_sink_create(const fc_perread_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts || !opts->output_root) return sink;

    perread_sink *state = calloc(1, sizeof(*state));
    state->include_alignment = opts->include_alignment;
    state->output_unmapped = opts->output_unmapped;
    char *path = fc_sink_join_path(opts->output_root, "per_read.tsv");
    ensure_parent_dir_exists(path);
    state->fp = fopen(path, "w");
    if (state->fp == NULL) {
        fprintf(stderr, "Error opening per-read file '%s' for writing.\n", path);
        exit(EXIT_FAILURE);
    }
    free(path);
    if (opts->sample && opts->sample[0]) {
        state->sample = strdup(opts->sample);
    }

    perread_sink_write_header(state);

    sink.ctx = state;
    sink.add_read = perread_sink_add_read;
    sink.destroy = perread_sink_destroy;
    return sink;
}
