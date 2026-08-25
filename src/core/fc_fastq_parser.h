#ifndef FASTCAT_FC_FASTQ_PARSER_H
#define FASTCAT_FC_FASTQ_PARSER_H

#include <stddef.h>

#include "htslib/sam.h"

#include "core/kseq_gz.h"

typedef struct fc_fastq_record {
    const char *name;
    size_t name_len;
    const char *comment;
    size_t comment_len;
    const char *sequence;
    const char *quality;
    size_t sequence_len;
} fc_fastq_record;

typedef struct fc_fastq_parser fc_fastq_parser;

fc_fastq_parser *fc_fastq_parser_create(fc_fastq_stream *stream);
void fc_fastq_parser_destroy(fc_fastq_parser *parser);

/* Return one for a record, zero at EOF, -2 for truncated quality, -4 for
 * missing quality, or -3 for malformed input/allocation/reader errors. The
 * returned field spans remain valid until the next call. */
int fc_fastq_parser_next(
    fc_fastq_parser *parser, bam1_t *bam, fc_fastq_record *record);

/* Build the unmapped BAM representation used for FASTQ input. Fields need not
 * be NUL-terminated; sequence and quality must each contain sequence_len
 * bytes. Qualities are ASCII FASTQ values rather than numeric Phred values. */
int fc_fastq_record_to_bam(bam1_t *bam,
    const char *name, size_t name_len,
    const char *comment, size_t comment_len,
    const char *sequence, const char *quality, size_t sequence_len);

#endif
