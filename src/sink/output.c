#include "fc_sinks.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "common.h"
#include "sink/common.h"
#include "version.h"

// default buffer size for writing to gzip, this is large enough that most reads will not require
// the write buffer to be resized. See gz_fprintf_bridge() below.
// The size is also used with gzbuffer() when opening gzFile handles.
#define GZBUFSIZE 131072  // 128 kB

typedef int (*sink_write_fn)(void *handle, const char *fmt, ...);

static int fprintf_bridge(void *h, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vfprintf((FILE *)h, fmt, ap);
    va_end(ap);
    return r;
}

static int gz_fprintf_bridge(void *h, const char *fmt, ...) {
    size_t bufsize = GZBUFSIZE;
    char *buf = (char *)xalloc(bufsize, sizeof(char), "gzbuffer");

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
        buf = (char *)xrealloc(buf, bufsize, "gzbuffer");
        va_start(ap, fmt);
        need = vsnprintf(buf, bufsize, fmt, ap);
        va_end(ap);
    }
    int w = gzwrite((gzFile)h, buf, (unsigned)need);
    free(buf);
    return (w == 0) ? -1 : need;
}

typedef struct {
    char *output;
    size_t reads_per_file;
    size_t barcode_count;
    bool reheader;
    bool write_bam;
    gzFile *handles;
    htsFile **bam_files;
    bam_hdr_t *bam_hdr;
    htsThreadPool hts_pool;
    size_t *reads_written;
    size_t *file_index;
} output_sink;

static void output_write_fastq(output_sink *sink, bam1_t *bam, read_meta meta, void *handle) {
    sink_write_fn sink_write = handle == stdout ? fprintf_bridge : gz_fprintf_bridge;

    static const char *fq_comment_fmt = "@%s %s\n%s\n+\n%s\n";
    static const char *sam_comment_fmt = "@%s\t%s\n%s\n+\n%s\n";
    static const char *no_comment_fmt = "@%s\n%s\n+\n%s\n";

    char *qname = bam_get_qname(bam);
    size_t len = bam->core.l_qseq;
    uint8_t *seq_data = bam_get_seq(bam);
    uint8_t *qual_data = bam_get_qual(bam);

    char *seq_str = xalloc(len + 1, sizeof(char), "bam seq");
    char *qual_str = xalloc(len + 1, sizeof(char), "bam qual");
    for (size_t i = 0; i < len; ++i) {
        seq_str[i] = seq_nt16_str[bam_seqi(seq_data, i)];
        qual_str[i] = (qual_data[0] == 0xff) ? '!' : (char)(qual_data[i] + 33);
    }
    seq_str[len] = '\0';
    qual_str[len] = '\0';

    if (sink->reheader) {
        if (meta->tags_str && meta->tags_str->l > 0) {
            sink_write(handle, sam_comment_fmt, qname, meta->tags_str->s, seq_str, qual_str);
        } else {
            sink_write(handle, no_comment_fmt, qname, seq_str, qual_str);
        }
    } else {
        const char *comment = (meta->orig_comment && meta->orig_comment[0]) ? meta->orig_comment : NULL;
        if (comment && comment[0]) {
            sink_write(handle, fq_comment_fmt, qname, comment, seq_str, qual_str);
        } else {
            sink_write(handle, no_comment_fmt, qname, seq_str, qual_str);
        }
    }

    free(seq_str);
    free(qual_str);
}

static void output_write_bam(output_sink *sink, bam1_t *bam, htsFile *handle) {
    if (bam == NULL) {
        fprintf(stderr, "Internal error: BAM record missing for output.\n");
        exit(1);
    }

    if (sam_write1(handle, sink->bam_hdr, bam) < 0) {
        fprintf(stderr, "Error writing read to BAM file.\n");
        exit(1);
    }
}

static void output_create_filepath(output_sink *sink, size_t barcode, char **path, char **filepath) {
    int ex_size = (sink->reads_per_file == 0) ? 0 : 5;
    char *extra = calloc(ex_size + 1, sizeof(char));
    if (sink->reads_per_file != 0) {
        sprintf(extra, "_%04zu", sink->file_index[barcode]);
    }

    if (barcode == 0) {
        *path = (char *)calloc(strlen(sink->output) + 15, sizeof(char));
        sprintf(*path, "%s/unclassified/", sink->output);
        *filepath = (char *)calloc(strlen(*path) + 22 + ex_size, sizeof(char));
        sprintf(*filepath, "%sunclassified%s.%s", *path, extra, sink->write_bam ? "bam" : "fastq.gz");
    } else {
        *path = (char *)calloc(strlen(sink->output) + 14, sizeof(char));
        sprintf(*path, "%s/barcode%04lu/", sink->output, barcode);
        *filepath = (char *)calloc(strlen(*path) + 21 + ex_size, sizeof(char));
        sprintf(*filepath, "%sbarcode%04lu%s.%s", *path, barcode, extra, sink->write_bam ? "bam" : "fastq.gz");
    }
    free(extra);
}

static void output_ensure_directory(output_sink *sink, size_t barcode, char *path) {
    if (sink->reads_per_file == 0
            || (sink->file_index[barcode] == 0 && sink->reads_written[barcode] == 0)) {
        int rtn = mkdir_hier(path);
        if (rtn == -1) {
            fprintf(stderr, "Failed to create barcode directory '%s\n'.", path);
            exit(1);
        }
    }
}

static void output_sink_add_read(void *ctx, const fc_read *read, const fc_read_ctx *read_ctx,
        const char *fname, const sam_hdr_t *hdr) {
    output_sink *sink = (output_sink *)ctx;
    (void)fname;
    (void)hdr;

    size_t barcode = read_ctx->barcode;
    if (barcode >= sink->barcode_count) {
        fprintf(stderr,
            "ERROR: Read's barcode number (%lu) exceeds configured barcode slots (%zu)\n",
            barcode, sink->barcode_count);
        exit(1);
    }

    if (sink->output == NULL) {
        if (sink->write_bam) {
            output_write_bam(sink, read->bam, sink->bam_files[0]);
        } else {
            output_write_fastq(sink, read->bam, read->meta, stdout);
        }
        return;
    }

    if (sink->reads_per_file != 0 && sink->reads_written[barcode] == sink->reads_per_file) {
        if (sink->write_bam) {
            hts_close(sink->bam_files[barcode]);
            sink->bam_files[barcode] = NULL;
        } else {
            if (sink->handles[barcode] == NULL) {
                fprintf(stderr, "Unexpected output file status encountered.");
                exit(1);
            }
            gzflush(sink->handles[barcode], Z_FINISH);
            gzclose(sink->handles[barcode]);
            sink->handles[barcode] = NULL;
        }
        sink->file_index[barcode]++;
        sink->reads_written[barcode] = 0;
    }

    {
        char *path = NULL;
        char *filepath = NULL;
        if (sink->write_bam) {
            if (sink->bam_files[barcode] == NULL) {
                output_create_filepath(sink, barcode, &path, &filepath);
                output_ensure_directory(sink, barcode, path);
                sink->bam_files[barcode] = hts_open(filepath, "wb");
                hts_set_opt(sink->bam_files[barcode], HTS_OPT_THREAD_POOL, &sink->hts_pool);
                if (sam_hdr_write(sink->bam_files[barcode], sink->bam_hdr)) {
                    fprintf(stderr, "Error writing header to BAM file\n");
                    exit(1);
                }
            }
            output_write_bam(sink, read->bam, sink->bam_files[barcode]);
        } else {
            if (sink->handles[barcode] == NULL) {
                output_create_filepath(sink, barcode, &path, &filepath);
                output_ensure_directory(sink, barcode, path);
                sink->handles[barcode] = gzopen(filepath, "wb");
                gzbuffer(sink->handles[barcode], GZBUFSIZE);
            }
            output_write_fastq(sink, read->bam, read->meta, sink->handles[barcode]);
        }
        free(filepath);
        free(path);
    }

    sink->reads_written[barcode]++;
}

static void output_sink_destroy(void *ctx) {
    output_sink *sink = (output_sink *)ctx;
    for (size_t i = 0; i < sink->barcode_count; ++i) {
        if (sink->write_bam) {
            if (sink->bam_files && sink->bam_files[i] != NULL) {
                hts_close(sink->bam_files[i]);
            }
        } else {
            if (sink->handles && sink->handles[i] != NULL) {
                gzflush(sink->handles[i], Z_FINISH);
                gzclose(sink->handles[i]);
            }
        }
    }
    if (sink->write_bam) {
        if (sink->hts_pool.pool) {
            hts_tpool_destroy(sink->hts_pool.pool);
        }
        free(sink->bam_files);
        bam_hdr_destroy(sink->bam_hdr);
    } else {
        free(sink->handles);
    }
    free(sink->reads_written);
    free(sink->file_index);
    free(sink->output);
    free(sink);
}

fc_sink fc_output_sink_create(const fc_output_sink_opts *opts) {
    fc_sink sink = {0};
    if (!opts) return sink;

    // Convention: non-NULL output_root means demuxed file output under <output_root>/demux; NULL streams to stdout.
    output_sink *state = calloc(1, sizeof(*state));
    state->output = NULL;
    if (opts->output_root) {
        char *path = fc_sink_join_path(opts->output_root, "demux");
        state->output = fc_sink_strip_path(path);
        free(path);
    }
    state->reads_per_file = opts->reads_per_file;
    state->barcode_count = opts->barcode_count ? opts->barcode_count : 1;
    state->reheader = opts->reheader;
    state->write_bam = opts->write_bam;
    state->reads_written = calloc(state->barcode_count, sizeof(size_t));
    state->file_index = calloc(state->barcode_count, sizeof(size_t));

    if (state->output != NULL) {
        int rtn = mkdir_hier(state->output);
        if (rtn == -1 && errno != EEXIST) {
            fprintf(stderr,
                "Error: Cannot create output directory '%s'. Check location is writeable and directory does not exist.\n",
                state->output);
            exit(EXIT_FAILURE);
        }
    }

    if (state->write_bam) {
        fprintf(stderr, "Using %d threads for BAM writing\n", opts->threads);
        state->hts_pool.pool = hts_tpool_init(opts->threads);
        state->hts_pool.qsize = 0;
        if (state->hts_pool.pool == NULL) {
            fprintf(stderr, "Error creating thread pool\n");
            exit(1);
        }
        state->bam_hdr = sam_hdr_init();
        sam_hdr_add_line(state->bam_hdr, "HD", "VN", SAM_FORMAT_VERSION, "SO", "unsorted", NULL);
        sam_hdr_add_line(state->bam_hdr, "PG", "ID", "fastcat", "PN", "fastcat", "VN", argp_program_version, NULL);
        state->bam_files = calloc(state->barcode_count, sizeof(htsFile *));
        if (state->output == NULL) {
            state->bam_files[0] = hts_open("-", "wb");
            hts_set_opt(state->bam_files[0], HTS_OPT_THREAD_POOL, &state->hts_pool);
            if (sam_hdr_write(state->bam_files[0], state->bam_hdr)) {
                fprintf(stderr, "Error writing header to BAM on stdout\n");
                exit(1);
            }
        }
    } else {
        state->handles = calloc(state->barcode_count, sizeof(gzFile));
    }

    sink.ctx = state;
    sink.add_read = output_sink_add_read;
    sink.destroy = output_sink_destroy;
    return sink;
}
