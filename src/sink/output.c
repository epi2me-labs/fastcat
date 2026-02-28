#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "../fc_log.h"
#include "common.h"
#include "fc_sinks.h"
#include "sink/common.h"
#include "version.h"

// default buffer size for writing to gzip, this is large enough that most reads will not require
// the write buffer to be resized. See gz_fprintf_bridge() below.
// The size is also used with gzbuffer() when opening gzFile handles.
#define GZBUFSIZE 131072 // 128 kB

typedef int (*sink_write_fn)(void *handle, const char *fmt, ...);

static int fprintf_bridge(void *h, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int r = vfprintf((FILE *)h, fmt, ap);
    va_end(ap);
    return r;
}

static int gz_fprintf_bridge(void *h, const char *fmt, ...)
{
    size_t bufsize = GZBUFSIZE;
    char *buf = fc_alloc(bufsize, sizeof(char), "gzbuffer");

    va_list ap;
    va_start(ap, fmt);
    int need = vsnprintf(buf, bufsize, fmt, ap);
    va_end(ap);

    if (need < 0) {
        free(buf);
        return -1;
    }

    if ((size_t)need >= bufsize) {
        bufsize = (size_t)need + 1;
        buf = fc_realloc(buf, bufsize, "gzbuffer");
        va_start(ap, fmt);
        need = vsnprintf(buf, bufsize, fmt, ap);
        va_end(ap);
    }
    int w = gzwrite((gzFile)h, buf, (unsigned)need);
    free(buf);
    return (w == 0) ? -1 : need;
}

typedef struct {
    fc_output_sink_opts cfg; // cached configuration
    char *output;
    gzFile *handles;
    htsFile **bam_files;
    bam_hdr_t *bam_hdr;
    htsThreadPool hts_pool;
    size_t *reads_written;
    size_t *file_index;
} output_sink;

static void output_write_fastq(output_sink *sink, const fc_read *read, void *handle)
{
    (void)sink;
    sink_write_fn sink_write = handle == stdout ? fprintf_bridge : gz_fprintf_bridge;

    static const char *sam_comment_fmt = "@%s\t%s\n%s\n+\n%s\n";
    static const char *no_comment_fmt = "@%s\n%s\n+\n%s\n";

    bam1_t *bam = read->bam;
    char *qname = bam_get_qname(bam);
    size_t len = bam->core.l_qseq;
    uint8_t *seq_data = bam_get_seq(bam);
    uint8_t *qual_data = bam_get_qual(bam);

    char *seq_str = fc_alloc(len + 1, sizeof(char), "bam seq");
    char *qual_str = fc_alloc(len + 1, sizeof(char), "bam qual");
    for (size_t i = 0; i < len; ++i) {
        seq_str[i] = seq_nt16_str[bam_seqi(seq_data, i)];
        qual_str[i] = (qual_data[0] == 0xff) ? '!' : (char)(qual_data[i] + 33);
    }
    seq_str[len] = '\0';
    qual_str[len] = '\0';

    int rc = 0;
    const kstring_t *tags = fc_read_get_tags_str(read);
    if (tags && tags->l > 0) {
        rc = sink_write(handle, sam_comment_fmt, qname, tags->s, seq_str, qual_str);
    } else {
        rc = sink_write(handle, no_comment_fmt, qname, seq_str, qual_str);
    }
    if (rc < 0) {
        FC_LOG_FATAL("Error writing FASTQ output.");
    }

    free(seq_str);
    free(qual_str);
}

static void output_write_bam(output_sink *sink, bam1_t *bam, htsFile *handle)
{
    if (bam == NULL) {
        FC_LOG_FATAL("Internal error: BAM record missing for output.");
    }

    if (sam_write1(handle, sink->bam_hdr, bam) < 0) {
        FC_LOG_FATAL("Error writing read to BAM file.");
    }
}

static void output_create_filepath(output_sink *sink, size_t barcode, char **path, char **filepath)
{
    char *extra = NULL;
    if (sink->cfg.reads_per_file != 0) {
        extra = fc_asprintf("output filename suffix", "_%04zu", sink->file_index[barcode]);
    } else {
        extra = fc_alloc(1, sizeof(char), "output filename suffix");
    }

    char stem_buf[32] = { 0 };
    const char *stem = fc_sink_barcode_stem(barcode, stem_buf, sizeof(stem_buf));
    if (!stem) {
        FC_LOG_FATAL("Failed to format barcode stem for output path.");
    }
    *path = fc_join_path(sink->output, stem);
    const char *ext = sink->cfg.write_bam ? "bam" : "fastq.gz";
    *filepath = fc_asprintf("output filepath", "%s/%s%s.%s", *path, stem, extra, ext);
    free(extra);
}

static void output_ensure_directory(output_sink *sink, size_t barcode, char *path)
{
    if (sink->cfg.reads_per_file == 0
        || (sink->file_index[barcode] == 0 && sink->reads_written[barcode] == 0)) {
        // Fail if barcode subdirectory already exists (fresh directory required)
        int rtn = mkdir_hier(path, true);
        if (rtn == -1) {
            FC_LOG_FATAL("Failed to create barcode directory '%s'.", path);
        }
    }
}

static void output_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
    const char *fname, const sam_hdr_t *hdr)
{
    output_sink *sink = (output_sink *)ctx;
    (void)fname;
    (void)hdr;

    size_t barcode = read_ctx->barcode;
    if (barcode >= sink->cfg.barcode_count) {
        FC_LOG_FATAL("Read's barcode number (%lu) exceeds configured barcode slots (%zu)",
            barcode, sink->cfg.barcode_count);
    }

    if (sink->output == NULL) {
        if (sink->cfg.write_bam) {
            output_write_bam(sink, read->bam, sink->bam_files[0]);
        } else {
            output_write_fastq(sink, read, stdout);
        }
        return;
    }

    if (sink->cfg.reads_per_file != 0 && sink->reads_written[barcode] == sink->cfg.reads_per_file) {
        if (sink->cfg.write_bam) {
            if (hts_close(sink->bam_files[barcode]) < 0) {
                FC_LOG_FATAL("Error closing BAM output file.");
            }
            sink->bam_files[barcode] = NULL;
        } else {
            if (sink->handles[barcode] == NULL) {
                FC_LOG_FATAL("Unexpected output file status encountered.");
            }
            if (gzclose(sink->handles[barcode]) != Z_OK) {
                FC_LOG_FATAL("Error closing FASTQ output file.");
            }
            sink->handles[barcode] = NULL;
        }
        sink->file_index[barcode]++;
        sink->reads_written[barcode] = 0;
    }

    {
        char *path = NULL;
        char *filepath = NULL;
        if (sink->cfg.write_bam) {
            if (sink->bam_files[barcode] == NULL) {
                output_create_filepath(sink, barcode, &path, &filepath);
                output_ensure_directory(sink, barcode, path);
                sink->bam_files[barcode] = hts_open(filepath, "wb");
                if (sink->bam_files[barcode] == NULL) {
                    FC_LOG_FATAL("Error opening BAM output file '%s'.", filepath);
                }
                if (hts_set_opt(sink->bam_files[barcode], HTS_OPT_THREAD_POOL, &sink->hts_pool) < 0) {
                    FC_LOG_FATAL("Error configuring BAM output file '%s'.", filepath);
                }
                if (sam_hdr_write(sink->bam_files[barcode], sink->bam_hdr)) {
                    FC_LOG_FATAL("Error writing header to BAM output file '%s'.", filepath);
                }
            }
            output_write_bam(sink, read->bam, sink->bam_files[barcode]);
        } else {
            if (sink->handles[barcode] == NULL) {
                output_create_filepath(sink, barcode, &path, &filepath);
                output_ensure_directory(sink, barcode, path);
                sink->handles[barcode] = gzopen(filepath, "wb");
                if (sink->handles[barcode] == NULL) {
                    FC_LOG_FATAL("Error opening FASTQ output file '%s'.", filepath);
                }
                if (gzbuffer(sink->handles[barcode], GZBUFSIZE) != 0) {
                    FC_LOG_FATAL("Error configuring FASTQ output file '%s'.", filepath);
                }
            }
            output_write_fastq(sink, read, sink->handles[barcode]);
        }
        free(filepath);
        free(path);
    }

    sink->reads_written[barcode]++;
}

static void output_sink_destroy(void *ctx)
{
    output_sink *sink = (output_sink *)ctx;
    for (size_t i = 0; i < sink->cfg.barcode_count; ++i) {
        if (sink->cfg.write_bam) {
            if (sink->bam_files && sink->bam_files[i] != NULL) {
                if (hts_close(sink->bam_files[i]) < 0) {
                    FC_LOG_FATAL("Error closing BAM output file.");
                }
            }
        } else {
            if (sink->handles && sink->handles[i] != NULL) {
                if (gzclose(sink->handles[i]) != Z_OK) {
                    FC_LOG_FATAL("Error closing FASTQ output file.");
                }
            }
        }
    }
    if (sink->cfg.write_bam) {
        if (sink->hts_pool.pool) {
            hts_tpool_destroy(sink->hts_pool.pool);
        }
        free(sink->bam_files);
        bam_hdr_destroy(sink->bam_hdr);
    } else {
        if (sink->output == NULL && fflush(stdout) == EOF) {
            FC_LOG_FATAL("Error flushing FASTQ stdout output.");
        }
        free(sink->handles);
    }
    free(sink->reads_written);
    free(sink->file_index);
    free(sink->output);
    free(sink);
}

fc_sink fc_output_sink_create(const fc_output_sink_opts *opts)
{
    fc_sink sink = { 0 };
    if (!opts) {
        return sink;
    }

    // Convention: non-NULL output_root means per-barcode output under <output_root>/{unclassified,barcodeXXXX}; NULL streams to stdout.
    output_sink *state = fc_alloc(1, sizeof(*state), "output sink");
    state->cfg = *opts;
    state->cfg.barcode_count = state->cfg.barcode_count ? state->cfg.barcode_count : 1;
    state->output = NULL;
    if (opts->output_root) {
        state->output = fc_strip_path(opts->output_root);
        if (!state->output) {
            FC_LOG_FATAL("Failed to prepare output path.");
        }
    }
    state->reads_written = calloc(state->cfg.barcode_count, sizeof(size_t));
    state->file_index = calloc(state->cfg.barcode_count, sizeof(size_t));
    if (!state->reads_written || !state->file_index) {
        FC_LOG_FATAL("Failed to allocate output sink counters.");
    }

    if (state->output != NULL) {
        // Allow directory to exist since it may have been created by fc_program already
        int rtn = mkdir_hier(state->output, false);
        if (rtn == -1) {
            FC_LOG_FATAL("Failed to create output directory '%s'.", state->output);
        }
    }

    if (state->cfg.write_bam) {
        FC_LOG_INFO("Using %d threads for BAM writing", opts->threads);
        state->hts_pool.pool = hts_tpool_init(opts->threads);
        state->hts_pool.qsize = 0;
        if (state->hts_pool.pool == NULL) {
            FC_LOG_FATAL("Error creating thread pool");
        }
        state->bam_hdr = sam_hdr_init();
        if (state->bam_hdr == NULL) {
            FC_LOG_FATAL("Error creating BAM header for output.");
        }
        sam_hdr_add_line(state->bam_hdr, "HD", "VN", SAM_FORMAT_VERSION, "SO", "unsorted", NULL);
        sam_hdr_add_line(state->bam_hdr, "PG", "ID", "fastcat", "PN", "fastcat", "VN", argp_program_version, NULL);
        state->bam_files = calloc(state->cfg.barcode_count, sizeof(htsFile *));
        if (state->bam_files == NULL) {
            FC_LOG_FATAL("Failed to allocate BAM output handles.");
        }
        if (state->output == NULL) {
            state->bam_files[0] = hts_open("-", "wb");
            if (state->bam_files[0] == NULL) {
                FC_LOG_FATAL("Error opening BAM stdout output.");
            }
            if (hts_set_opt(state->bam_files[0], HTS_OPT_THREAD_POOL, &state->hts_pool) < 0) {
                FC_LOG_FATAL("Error configuring BAM stdout output.");
            }
            if (sam_hdr_write(state->bam_files[0], state->bam_hdr)) {
                FC_LOG_FATAL("Error writing header to BAM on stdout");
            }
        }
    } else {
        state->handles = calloc(state->cfg.barcode_count, sizeof(gzFile));
        if (state->handles == NULL) {
            FC_LOG_FATAL("Failed to allocate FASTQ output handles.");
        }
    }

    sink.ctx = state;
    sink.add_read = output_sink_add_read;
    sink.destroy = output_sink_destroy;
    return sink;
}
