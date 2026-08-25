#ifndef FASTCAT_KSEQ_GZ_H
#define FASTCAT_KSEQ_GZ_H

#include <stdbool.h>
#include <stdint.h>

#include "htslib/kseq.h"

typedef struct fc_fastq_stream fc_fastq_stream;

typedef struct {
    const unsigned char *data;
    size_t size;
} fc_fastq_span;

fc_fastq_stream *fc_fastq_stream_open_path(const char *path, unsigned workers);
fc_fastq_stream *fc_fastq_stream_open_fd(int fd);
int fc_fastq_stream_read(fc_fastq_stream *stream, void *buffer, int capacity);
int fc_fastq_stream_close(fc_fastq_stream *stream);
bool fc_fastq_stream_input_offset(const fc_fastq_stream *stream, uint64_t *out);
int fc_fastq_stream_next(fc_fastq_stream *stream, fc_fastq_span *span);

#ifndef KSEQ_DECLARED
KSEQ_DECLARE(fc_fastq_stream *)
#define KSEQ_DECLARED
#endif

#endif
