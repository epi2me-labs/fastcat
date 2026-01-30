#include "fc_sinks.h"

#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "fc_filter.h"
#include "sink/common.h"
#include "../fc_sinks.h"
#include "../fc_log.h"


#define FLAG_COLS 9

typedef struct flag_stats {
    size_t n_refs;
    size_t **counts;
    size_t *unmapped;
} flag_stats;

// Allocate and destroy flag statistics structures.
static flag_stats* create_flag_stats(size_t n_refs, bool store_unmapped) {
    flag_stats* stats = xalloc(1, sizeof(flag_stats), "flagstat");
    stats->n_refs = n_refs;
    stats->counts = xalloc(n_refs, sizeof(size_t*), "flagstat");
    stats->unmapped = store_unmapped ? xalloc(FLAG_COLS, sizeof(size_t), "flagstat") : NULL;

    for (size_t i = 0; i < n_refs; i++) {
        stats->counts[i] = xalloc(FLAG_COLS, sizeof(size_t), "flagstat");
    }
    return stats;
}

static void destroy_flag_stats(flag_stats* stats) {
    if (!stats) return;
    for (size_t i = 0; i < stats->n_refs; i++) {
        free(stats->counts[i]);
    }
    free(stats->counts);
    free(stats->unmapped);
    free(stats);
}

// counting alignment flags: total, primary, secondary, supplementary, unmapped, qcfail, duplicate, duplex, duplex_forming
static const size_t FLAG_MASK[8] = {
    0, 0, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP, BAM_FQCFAIL, BAM_FDUP, 0
};

static inline void flagstats_add_counts(const bam1_t* b, size_t* counts, const int duplex_code) {
    counts[0] += 1;
    counts[1] += fc_is_primary(b);
    for (size_t i = 2; i < 6; ++i) {
        counts[i] += ((b->core.flag & FLAG_MASK[i]) != 0);
    }
    counts[7] += (duplex_code == 1);
    counts[8] += (duplex_code == -1);
}

typedef struct {
    flag_stats* stats;
    bool region_mode;
    bool track_unmapped;
    const bam_hdr_t *hdr;
    FILE *fp;
    char *sample;
} flagstats_sink;

static void flagstats_sink_add_read(void* ctx, const fc_read* fr, const fc_read_ctx* read_ctx,
        const char* fname, const sam_hdr_t* hdr) {
    flagstats_sink* sink = (flagstats_sink*)ctx;
    (void)fname;
    (void)hdr;
    (void)read_ctx;
    if (!sink->stats) return;

    const bam1_t* b = fr->bam;
    bool is_unmapped = fc_is_unmapped(b);
    if (!is_unmapped) {
        // Protect against invalid tid values (e.g., malformed NM test cases).
        size_t tid = (size_t)b->core.tid;
        if (!sink->region_mode && (tid >= sink->stats->n_refs)) {
            return;
        }
        size_t* counts = sink->region_mode ? sink->stats->counts[0] : sink->stats->counts[tid];
        flagstats_add_counts(b, counts, fr->meta->duplex_code);
    }

    if (sink->track_unmapped && sink->stats->unmapped && is_unmapped) {
        flagstats_add_counts(b, sink->stats->unmapped, fr->meta->duplex_code);
    }
}

static void flagstats_write_header(FILE* fh, const char* sample) {
    if (!fh) return;
    const char *sample_col = sample ? "\tsample_name" : "";
    fprintf(fh,
        "ref%s\ttotal\tprimary\tsecondary\tsupplementary\tunmapped\tqcfail\tduplicate\tduplex\tduplex_forming\n",
        sample_col);
}

static void flagstats_write_line(FILE* fh, const char* sample, const char* label, size_t* stats) {
    if (!fh) return;
    if (sample) {
        fprintf(fh,
            "%s\t%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
            label, sample,
            stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7], stats[8]);
    } else {
        fprintf(fh,
            "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
            label,
            stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7], stats[8]);
    }
}

static void flagstats_sink_finish(void* ctx) {
    flagstats_sink* sink = (flagstats_sink*)ctx;
    if (!sink || !sink->stats || !sink->fp) return;

    for (size_t i = 0; i < sink->stats->n_refs; ++i) {
        const char *label = NULL;
        if (sink->hdr && i < (size_t)sink->hdr->n_targets) {
            label = sam_hdr_tid2name((sam_hdr_t *)sink->hdr, (int)i);
        }
        char tmp[32];
        if (!label) {
            snprintf(tmp, sizeof(tmp), "%zu", i);
            label = tmp;
        }
        flagstats_write_line(sink->fp, sink->sample, label, sink->stats->counts[i]);
    }
    if (sink->track_unmapped && sink->stats->unmapped) {
        flagstats_write_line(sink->fp, sink->sample, "*", sink->stats->unmapped);
    }
}

static void flagstats_sink_destroy(void *ctx) {
    flagstats_sink *sink = (flagstats_sink *)ctx;
    if (sink->fp) fclose(sink->fp);
    free(sink->sample);
    destroy_flag_stats(sink->stats);
    free(sink);
}

fc_sink fc_flagstats_sink_create(const fc_flagstats_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts || !opts->output_root) return sink;

    flagstats_sink* state = calloc(1, sizeof(*state));
    size_t n_refs = opts->region_mode ? 1 : (opts->hdr ? (size_t)opts->hdr->n_targets : 0);
    state->stats = create_flag_stats(n_refs, opts->track_unmapped);
    state->region_mode = opts->region_mode;
    state->track_unmapped = opts->track_unmapped;
    state->hdr = opts->hdr;
    if (opts->sample && opts->sample[0]) {
        state->sample = strdup(opts->sample);
    }
    if (opts->output_root) {
        char *path = fc_sink_join_path(opts->output_root, "flagstats.tsv");
        ensure_parent_dir_exists(path);
        state->fp = fopen(path, "w");
        if (!state->fp) {
            FC_LOG_ERROR("Cannot open file '%s' for writing.", path);
            exit(EXIT_FAILURE);
        }
        flagstats_write_header(state->fp, state->sample);
        free(path);
    }

    sink.ctx = state;
    sink.add_read = flagstats_sink_add_read;
    sink.finish = flagstats_sink_finish;
    sink.destroy = flagstats_sink_destroy;
    return sink;
}
