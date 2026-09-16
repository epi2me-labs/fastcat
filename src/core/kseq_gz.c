#define KSEQ_DECLARED
#include "core/kseq_gz.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>

#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
#include "zngp.h"
#endif

struct fc_fastq_stream {
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    ZNGP *reader;
#else
    gzFile gz;
#endif
};

#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
enum { FC_HIGH_THREAD_VIEW_BUFFER_BYTES = 32 * 1024 * 1024 };
#endif

fc_fastq_stream *fc_fastq_stream_open_path(const char *path, unsigned workers)
{
    if (!path) {
        errno = EINVAL;
        return NULL;
    }
    fc_fastq_stream *stream = calloc(1, sizeof(*stream));
    if (!stream) {
        return NULL;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    stream->reader = zngp_open(path, "rb");
    if (!stream->reader) {
        free(stream);
        return NULL;
    }
    if (zngp_mt(stream->reader, (int)(workers ? workers : 1), 0) != 0) {
        int saved_errno = errno;
        (void)zngp_close(stream->reader);
        free(stream);
        errno = saved_errno;
        return NULL;
    }
    if (workers > 16 && zngp_set_view_buffer(stream->reader, FC_HIGH_THREAD_VIEW_BUFFER_BYTES) != 0) {
        int saved_errno = errno;
        (void)zngp_close(stream->reader);
        free(stream);
        errno = saved_errno;
        return NULL;
    }
    return stream;
#else
    (void)workers;
    stream->gz = gzopen(path, "r");
    if (!stream->gz) {
        free(stream);
        return NULL;
    }
    return stream;
#endif
}

fc_fastq_stream *fc_fastq_stream_open_fd(int fd, unsigned workers)
{
    fc_fastq_stream *stream = calloc(1, sizeof(*stream));
    if (!stream) {
        return NULL;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    stream->reader = zngp_dopen(fd, "rb");
    if (!stream->reader) {
        free(stream);
        return NULL;
    }
    if (zngp_mt(stream->reader, (int)(workers ? workers : 1), 0) != 0) {
        int saved_errno = errno;
        (void)zngp_close(stream->reader);
        free(stream);
        errno = saved_errno;
        return NULL;
    }
#else
    (void)workers;
    stream->gz = gzdopen(fd, "r");
    if (!stream->gz) {
        free(stream);
        return NULL;
    }
#endif
    return stream;
}

int fc_fastq_stream_read(fc_fastq_stream *stream, void *buffer, int capacity)
{
    if (!stream) {
        return -1;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    return (int)zngp_read(stream->reader, buffer, (size_t)capacity);
#else
    return gzread(stream->gz, buffer, (unsigned)capacity);
#endif
}

int fc_fastq_stream_close(fc_fastq_stream *stream)
{
    if (!stream) {
        return -1;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    int status = zngp_close(stream->reader);
#else
    int status = gzclose(stream->gz) == Z_OK ? 0 : -1;
#endif
    free(stream);
    return status;
}

bool fc_fastq_stream_input_offset(const fc_fastq_stream *stream,
    uint64_t *out)
{
    if (!stream || !out) {
        return false;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    if (!stream->reader) {
        return false;
    }
    off_t offset = zngp_offset(stream->reader);
#else
    if (!stream->gz) {
        return false;
    }
    z_off_t offset = gzoffset(stream->gz);
#endif
    if (offset < 0) {
        return false;
    }
    *out = (uint64_t)offset;
    return true;
}

int fc_fastq_stream_next(fc_fastq_stream *stream, fc_fastq_span *span)
{
    if (!stream || !span) {
        errno = EINVAL;
        return -1;
    }
#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    zngp_view view;
    int status = zngp_next(stream->reader, &view);
    span->data = view.data;
    span->size = view.size;
    return status;
#else
    span->data = NULL;
    span->size = 0;
    errno = ENOTSUP;
    return -1;
#endif
}

// Match purrflate's reader buffer to avoid waking its producer for every
// 16 KiB kseq refill. This affects only Fastcat's kseq instantiation.
enum { FC_KSEQ_BUFFER_BYTES = 1024 * 1024 };
KSTREAM_INIT2(, fc_fastq_stream *, fc_fastq_stream_read, FC_KSEQ_BUFFER_BYTES)
__KSEQ_TYPE(fc_fastq_stream *)
__KSEQ_BASIC(, fc_fastq_stream *)
__KSEQ_READ()
