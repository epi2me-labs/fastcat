#include "core/fc_fastq_parser.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__i386__) || defined(__x86_64__)
#include <immintrin.h>
#endif

#include "htslib/hts.h"

#include "core/fc_fastqcomments.h"

#if defined(__i386__) || defined(__x86_64__)
// VBMI provides a full-width byte lookup, avoiding AVX2's independent
// 128-bit shuffle lanes and the corresponding result compaction.
__attribute__((target("avx512f,avx512bw,avx512vbmi"))) static size_t pack_fastq_sequence_avx512vbmi(
    uint8_t *packed, const uint8_t *ascii, size_t length)
{
    const __m512i zero = _mm512_setzero_si512();
    // Index by the low six ASCII bits. Separate expected-byte lookup retains
    // exact validation for uppercase and lowercase A/C/G/T/N.
    static const uint8_t code_lut_bytes[64] = {
        0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0,
        0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0,
        0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    static const uint8_t expected_lut_bytes[64] = {
        0, 'A', 0, 'C', 0, 0, 0, 'G', 0, 0, 0, 0, 0, 0, 'N', 0,
        0, 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 'a', 0, 'c', 0, 0, 0, 'g', 0, 0, 0, 0, 0, 0, 'n', 0,
        0, 0, 0, 0, 't', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    static const uint8_t even_indices[64] = {
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
        32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    static const uint8_t odd_indices[64] = {
        1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
        33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    const __m512i code_lut = _mm512_loadu_si512(code_lut_bytes);
    const __m512i expected_lut = _mm512_loadu_si512(expected_lut_bytes);
    const __m512i even_bytes = _mm512_loadu_si512(even_indices);
    const __m512i odd_bytes = _mm512_loadu_si512(odd_indices);

    size_t i = 0;
    for (; i + 64 <= length; i += 64) {
        __m512i bases = _mm512_loadu_si512((const void *)(ascii + i));
        __m512i codes = _mm512_permutexvar_epi8(bases, code_lut);
        __m512i expected = _mm512_permutexvar_epi8(bases, expected_lut);
        __mmask64 expected_matches = _mm512_cmpeq_epi8_mask(bases, expected);
        __mmask64 code_is_nonzero = _mm512_cmpneq_epi8_mask(codes, zero);

        if ((uint64_t)(expected_matches & code_is_nonzero) != UINT64_MAX) {
            for (size_t j = 0; j < 64; j += 2) {
                packed[(i + j) / 2] = (uint8_t)((seq_nt16_table[ascii[i + j]] << 4)
                    | seq_nt16_table[ascii[i + j + 1]]);
            }
            continue;
        }

        __m512i even = _mm512_permutexvar_epi8(even_bytes, codes);
        __m512i odd = _mm512_permutexvar_epi8(odd_bytes, codes);
        __m512i pairs = _mm512_or_si512(_mm512_slli_epi16(even, 4), odd);
        _mm256_storeu_si256(
            (__m256i *)(packed + i / 2), _mm512_castsi512_si256(pairs));
    }
    return i;
}

// AVX2 uses the same nibble lookup as the SSSE3 path, duplicated in each
// 128-bit lane. The packed halves are joined before the contiguous store.
__attribute__((target("avx2"))) static size_t pack_fastq_sequence_avx2(uint8_t *packed, const uint8_t *ascii, size_t length)
{
    const __m256i uppercase_mask = _mm256_set1_epi8((char)0xdf);
    const __m256i low_nibble_mask = _mm256_set1_epi8(0x0f);
    const __m256i high_nibble_mask = _mm256_set1_epi8((char)0xf0);
    const __m256i zero = _mm256_setzero_si256();
    const __m256i code_lut = _mm256_setr_epi8(
        0, 1, 0, 2, 8, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0,
        0, 1, 0, 2, 8, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0);
    const __m256i high_nibble_lut = _mm256_setr_epi8(
        0, 0x40, 0, 0x40, 0x50, 0, 0, 0x40,
        0, 0, 0, 0, 0, 0, 0x40, 0,
        0, 0x40, 0, 0x40, 0x50, 0, 0, 0x40,
        0, 0, 0, 0, 0, 0, 0x40, 0);
    const __m256i even_bytes = _mm256_setr_epi8(
        0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i odd_bytes = _mm256_setr_epi8(
        1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1,
        1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);

    size_t i = 0;
    for (; i + 32 <= length; i += 32) {
        __m256i bases = _mm256_loadu_si256((const __m256i *)(ascii + i));
        __m256i upper = _mm256_and_si256(bases, uppercase_mask);
        __m256i low_nibbles = _mm256_and_si256(upper, low_nibble_mask);
        __m256i codes = _mm256_shuffle_epi8(code_lut, low_nibbles);
        __m256i expected_high = _mm256_shuffle_epi8(high_nibble_lut, low_nibbles);
        __m256i high_matches = _mm256_cmpeq_epi8(
            _mm256_and_si256(upper, high_nibble_mask), expected_high);
        __m256i code_is_zero = _mm256_cmpeq_epi8(codes, zero);
        __m256i valid = _mm256_andnot_si256(code_is_zero, high_matches);

        if ((uint32_t)_mm256_movemask_epi8(valid) != UINT32_MAX) {
            for (size_t j = 0; j < 32; j += 2) {
                packed[(i + j) / 2] = (uint8_t)((seq_nt16_table[ascii[i + j]] << 4)
                    | seq_nt16_table[ascii[i + j + 1]]);
            }
            continue;
        }

        __m256i even = _mm256_shuffle_epi8(codes, even_bytes);
        __m256i odd = _mm256_shuffle_epi8(codes, odd_bytes);
        __m256i pairs = _mm256_or_si256(_mm256_slli_epi16(even, 4), odd);
        __m128i lower = _mm256_castsi256_si128(pairs);
        __m128i upper_pairs = _mm256_extracti128_si256(pairs, 1);
        __m128i contiguous = _mm_unpacklo_epi64(lower, upper_pairs);
        _mm_storeu_si128((__m128i *)(packed + i / 2), contiguous);
    }
    return i;
}

// Most basecalled reads contain only A/C/G/T/N. Pack those bases 16 at a time,
// while retaining the scalar table fallback for the complete BAM alphabet.
__attribute__((target("ssse3"))) static size_t pack_fastq_sequence_ssse3(uint8_t *packed, const uint8_t *ascii, size_t length)
{
    const __m128i uppercase_mask = _mm_set1_epi8((char)0xdf);
    const __m128i low_nibble_mask = _mm_set1_epi8(0x0f);
    const __m128i high_nibble_mask = _mm_set1_epi8((char)0xf0);
    const __m128i zero = _mm_setzero_si128();
    const __m128i code_lut = _mm_setr_epi8(
        0, 1, 0, 2, 8, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0);
    const __m128i high_nibble_lut = _mm_setr_epi8(
        0, 0x40, 0, 0x40, 0x50, 0, 0, 0x40,
        0, 0, 0, 0, 0, 0, 0x40, 0);
    const __m128i even_bytes = _mm_setr_epi8(
        0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i odd_bytes = _mm_setr_epi8(
        1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);

    size_t i = 0;
    for (; i + 16 <= length; i += 16) {
        __m128i bases = _mm_loadu_si128((const __m128i *)(ascii + i));
        __m128i upper = _mm_and_si128(bases, uppercase_mask);
        __m128i low_nibbles = _mm_and_si128(upper, low_nibble_mask);
        __m128i codes = _mm_shuffle_epi8(code_lut, low_nibbles);
        __m128i expected_high = _mm_shuffle_epi8(high_nibble_lut, low_nibbles);
        __m128i high_matches = _mm_cmpeq_epi8(
            _mm_and_si128(upper, high_nibble_mask), expected_high);
        __m128i code_is_zero = _mm_cmpeq_epi8(codes, zero);
        __m128i valid = _mm_andnot_si128(code_is_zero, high_matches);

        if (_mm_movemask_epi8(valid) != 0xffff) {
            for (size_t j = 0; j < 16; j += 2) {
                packed[(i + j) / 2] = (uint8_t)((seq_nt16_table[ascii[i + j]] << 4)
                    | seq_nt16_table[ascii[i + j + 1]]);
            }
            continue;
        }

        __m128i even = _mm_shuffle_epi8(codes, even_bytes);
        __m128i odd = _mm_shuffle_epi8(codes, odd_bytes);
        __m128i pairs = _mm_or_si128(_mm_slli_epi16(even, 4), odd);
        _mm_storel_epi64((__m128i *)(packed + i / 2), pairs);
    }
    return i;
}

#endif

static void pack_fastq_sequence(uint8_t *packed, const uint8_t *ascii, size_t length)
{
    size_t i = 0;
#if defined(__i386__) || defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f")
        && __builtin_cpu_supports("avx512bw")
        && __builtin_cpu_supports("avx512vbmi")) {
        i = pack_fastq_sequence_avx512vbmi(packed, ascii, length);
    }
    if (__builtin_cpu_supports("avx2")) {
        i += pack_fastq_sequence_avx2(
            packed + i / 2, ascii + i, length - i);
    }
    if (__builtin_cpu_supports("ssse3")) {
        i += pack_fastq_sequence_ssse3(
            packed + i / 2, ascii + i, length - i);
    }
#endif
    for (; i + 1 < length; i += 2) {
        packed[i / 2] = (uint8_t)((seq_nt16_table[ascii[i]] << 4) | seq_nt16_table[ascii[i + 1]]);
    }
    if (i < length) {
        packed[i / 2] = (uint8_t)(seq_nt16_table[ascii[i]] << 4);
    }
}

static void copy_fastq_qualities(
    uint8_t *restrict bam_qual, const uint8_t *restrict ascii_qual, size_t length)
{
    for (size_t i = 0; i < length; ++i) {
        bam_qual[i] = (uint8_t)(ascii_qual[i] - 33);
    }
}

// Constructor specialised for Fastcat's unmapped FASTQ input. Qualities are
// FASTQ ASCII rather than numeric Phred values.
static inline int fastcat_ubam_set1(
    bam1_t *bam,
    size_t name_len, const char *qname,
    size_t l_seq, const char *seq, const char *qual)
{
    if (!bam || !qname || !seq || !qual || name_len == 0 || name_len >= 254
        || l_seq > INT32_MAX) {
        errno = EINVAL;
        return -1;
    }

    size_t l_qname = name_len + 1;
    size_t qname_padding = 4 - l_qname % 4;
    size_t packed_seq_bytes = (l_seq + 1) / 2;
    size_t data_len = l_qname + qname_padding + packed_seq_bytes + l_seq;
    if (data_len > INT32_MAX
        || hts_resize(uint8_t, data_len, &bam->m_data, &bam->data, 0) < 0) {
        return -1;
    }

    memset(&bam->core, 0, sizeof(bam->core));
    bam->core.tid = -1;
    bam->core.pos = -1;
    bam->core.bin = hts_reg2bin(-1, 0, 14, 5);
    bam->core.flag = BAM_FUNMAP;
    bam->core.l_qname = (uint16_t)(l_qname + qname_padding);
    bam->core.l_extranul = (uint8_t)(qname_padding - 1);
    bam->core.l_qseq = (int32_t)l_seq;
    bam->core.mtid = -1;
    bam->core.mpos = -1;
    bam->l_data = (int)data_len;

    memcpy(bam->data, qname, name_len);
    memset(bam->data + name_len, 0, qname_padding + 1);

    uint8_t *packed_seq = bam_get_seq(bam);
    pack_fastq_sequence(packed_seq, (const uint8_t *)seq, l_seq);

    uint8_t *bam_qual = bam_get_qual(bam);
    copy_fastq_qualities(bam_qual, (const uint8_t *)qual, l_seq);

    return (int)data_len;
}

int fc_fastq_record_to_bam(bam1_t *bam,
    const char *name, size_t name_len,
    const char *comment, size_t comment_len,
    const char *seq, const char *qual, size_t seq_len)
{
    if (!bam || !name || !seq || !qual
        || (comment_len != 0 && !comment)) {
        errno = EINVAL;
        return -1;
    }
    if (fastcat_ubam_set1(bam, name_len, name, seq_len, seq, qual) < 0) {
        return -1;
    }

    kstring_t comment_view = {
        .l = comment_len,
        .m = comment_len,
        .s = (char *)comment
    };
    return comment_len == 0 || fastq_comment_to_bam(&comment_view, bam) == 0 ? 0 : -1;
}

typedef enum {
    PARSE_COMPLETE,
    PARSE_NEED_MORE,
    PARSE_QUAL_MISSING,
    PARSE_QUAL_ERROR,
    PARSE_INVALID
} parse_status;

typedef struct {
    size_t consumed;
    size_t header_end;
    size_t sequence_start;
    size_t plus_start;
    size_t quality_start;
    size_t sequence_len;
    size_t sequence_lines;
    size_t quality_lines;
} record_layout;

struct fc_fastq_parser {
    fc_fastq_stream *stream;
    fc_fastq_span span;
    size_t position;
    int span_active;
    int eof;
    char *carry;
    size_t carry_len;
    size_t carry_cap;
    char *sequence;
    size_t sequence_cap;
    char *quality;
    size_t quality_cap;
};

static int reserve(char **buffer, size_t *capacity, size_t required)
{
    if (required <= *capacity) {
        return 0;
    }
    size_t grown = *capacity ? *capacity : 4096;
    while (grown < required) {
        if (grown > SIZE_MAX / 2) {
            grown = required;
            break;
        }
        grown *= 2;
    }
    char *next = realloc(*buffer, grown);
    if (!next) {
        return -1;
    }
    *buffer = next;
    *capacity = grown;
    return 0;
}

static int append_carry(fc_fastq_parser *reader, const char *data, size_t size)
{
    if (size > SIZE_MAX - reader->carry_len
        || reserve(&reader->carry, &reader->carry_cap, reader->carry_len + size) != 0) {
        return -1;
    }
    memcpy(reader->carry + reader->carry_len, data, size);
    reader->carry_len += size;
    return 0;
}

static int next_line(const char *data, size_t size, size_t offset, int eof,
    size_t *line_len, size_t *next)
{
    if (offset >= size) {
        return 0;
    }
    const char *newline = memchr(data + offset, '\n', size - offset);
    if (!newline && !eof) {
        return 0;
    }
    size_t end = newline ? (size_t)(newline - data) : size;
    *line_len = end - offset;
    if (*line_len != 0 && data[end - 1] == '\r') {
        (*line_len)--;
    }
    *next = newline ? end + 1 : size;
    return 1;
}

static const char *find_header(const char *data, size_t size)
{
    if (size != 0 && (data[0] == '@' || data[0] == '>')) {
        return data;
    }
    const char *fastq = memchr(data, '@', size);
    const char *fasta = memchr(data, '>', size);
    if (!fastq) {
        return fasta;
    }
    return fasta && fasta < fastq ? fasta : fastq;
}

static parse_status parse_layout(
    const char *data, size_t size, int eof, record_layout *layout)
{
    memset(layout, 0, sizeof(*layout));
    size_t line_len;
    size_t next;
    if (!next_line(data, size, 0, eof, &line_len, &next)) {
        return PARSE_NEED_MORE;
    }
    if (line_len == 0 || (data[0] != '@' && data[0] != '>')) {
        return PARSE_INVALID;
    }
    layout->header_end = line_len;
    layout->sequence_start = next;

    size_t offset = next;
    for (;;) {
        size_t line_start = offset;
        if (!next_line(data, size, offset, eof, &line_len, &next)) {
            if (!eof) {
                return PARSE_NEED_MORE;
            }
            layout->consumed = size;
            return PARSE_QUAL_MISSING;
        }
        if (line_len != 0
            && (data[line_start] == '@' || data[line_start] == '>')) {
            layout->consumed = line_start;
            return PARSE_QUAL_MISSING;
        }
        if (line_len != 0 && data[line_start] == '+') {
            layout->plus_start = line_start;
            layout->quality_start = next;
            break;
        }
        if (line_len > SIZE_MAX - layout->sequence_len) {
            return PARSE_INVALID;
        }
        layout->sequence_len += line_len;
        layout->sequence_lines++;
        offset = next;
    }
    offset = layout->quality_start;
    size_t quality_len = 0;
    do {
        if (!next_line(data, size, offset, eof, &line_len, &next)) {
            if (!eof) {
                return PARSE_NEED_MORE;
            }
            layout->consumed = size;
            return PARSE_QUAL_ERROR;
        }
        if (line_len > SIZE_MAX - quality_len) {
            return PARSE_INVALID;
        }
        quality_len += line_len;
        layout->quality_lines++;
        offset = next;
        layout->consumed = offset;
        if (quality_len > layout->sequence_len) {
            return PARSE_QUAL_ERROR;
        }
    } while (quality_len < layout->sequence_len);

    if (quality_len == 0) {
        return PARSE_QUAL_MISSING;
    }
    layout->consumed = offset;
    return PARSE_COMPLETE;
}

static int copy_lines(char *output, const char *data,
    size_t start, size_t end, size_t *written)
{
    size_t offset = start;
    *written = 0;
    while (offset < end) {
        const char *newline = memchr(data + offset, '\n', end - offset);
        size_t line_end = newline ? (size_t)(newline - data) : end;
        size_t length = line_end - offset;
        if (length != 0 && data[line_end - 1] == '\r') {
            length--;
        }
        memcpy(output + *written, data + offset, length);
        *written += length;
        offset = newline ? line_end + 1 : end;
    }
    return 0;
}

static int build_record(fc_fastq_parser *reader, bam1_t *bam,
    fc_fastq_record *record, const char *data, const record_layout *layout)
{
    const char *header = data + 1;
    const char *header_end = data + layout->header_end;
    const char *name_end = header;
    while (name_end < header_end && *name_end != ' ' && *name_end != '\t') {
        name_end++;
    }
    const char *comment = name_end;
    while (comment < header_end && (*comment == ' ' || *comment == '\t')) {
        comment++;
    }

    const char *sequence = data + layout->sequence_start;
    const char *quality = data + layout->quality_start;
    if (layout->sequence_lines != 1 || layout->quality_lines != 1) {
        if (reserve(&reader->sequence, &reader->sequence_cap,
                layout->sequence_len + 1)
                != 0
            || reserve(&reader->quality, &reader->quality_cap,
                   layout->sequence_len + 1)
                != 0) {
            return -1;
        }
        size_t sequence_written;
        size_t quality_written;
        copy_lines(reader->sequence, data,
            layout->sequence_start, layout->plus_start, &sequence_written);
        copy_lines(reader->quality, data,
            layout->quality_start, layout->consumed, &quality_written);
        if (sequence_written != layout->sequence_len
            || quality_written != layout->sequence_len) {
            return -1;
        }
        reader->sequence[layout->sequence_len] = '\0';
        reader->quality[layout->sequence_len] = '\0';
        sequence = reader->sequence;
        quality = reader->quality;
    }

    record->name = header;
    record->name_len = (size_t)(name_end - header);
    record->comment = comment;
    record->comment_len = (size_t)(header_end - comment);
    record->sequence = sequence;
    record->quality = quality;
    record->sequence_len = layout->sequence_len;
    return fc_fastq_record_to_bam(bam,
        record->name, record->name_len,
        record->comment, record->comment_len,
        record->sequence, record->quality, record->sequence_len);
}

static int acquire_span(fc_fastq_parser *reader)
{
    if (reader->span_active && reader->position < reader->span.size) {
        return 1;
    }
    reader->span = (fc_fastq_span) { 0 };
    reader->position = 0;
    reader->span_active = 0;
    int status = fc_fastq_stream_next(reader->stream, &reader->span);
    if (status == 1) {
        reader->span_active = 1;
    } else if (status == 0) {
        reader->eof = 1;
    }
    return status;
}

static void consume_carry(fc_fastq_parser *reader, size_t consumed)
{
    if (consumed >= reader->carry_len) {
        reader->carry_len = 0;
        return;
    }
    memmove(reader->carry, reader->carry + consumed, reader->carry_len - consumed);
    reader->carry_len -= consumed;
}

fc_fastq_parser *fc_fastq_parser_create(fc_fastq_stream *stream)
{
    if (!stream) {
        errno = EINVAL;
        return NULL;
    }
    fc_fastq_parser *reader = calloc(1, sizeof(*reader));
    if (reader) {
        reader->stream = stream;
    }
    return reader;
}

void fc_fastq_parser_destroy(fc_fastq_parser *reader)
{
    if (!reader) {
        return;
    }
    free(reader->carry);
    free(reader->sequence);
    free(reader->quality);
    free(reader);
}

int fc_fastq_parser_next(
    fc_fastq_parser *reader, bam1_t *bam, fc_fastq_record *record)
{
    if (!reader || !bam || !record) {
        errno = EINVAL;
        return -3;
    }
    memset(record, 0, sizeof(*record));
    for (;;) {
        if (reader->carry_len != 0) {
            record_layout layout;
            parse_status parsed = parse_layout(
                reader->carry, reader->carry_len, reader->eof, &layout);
            if (parsed == PARSE_COMPLETE) {
                if (layout.consumed != reader->carry_len
                    || build_record(reader, bam, record, reader->carry, &layout) != 0) {
                    return -3;
                }
                reader->carry_len = 0;
                return 1;
            }
            if (parsed == PARSE_QUAL_MISSING || parsed == PARSE_QUAL_ERROR) {
                consume_carry(reader, layout.consumed);
                return parsed == PARSE_QUAL_MISSING ? -4 : -2;
            }
            if (parsed == PARSE_INVALID) {
                return -3;
            }
            int acquired = acquire_span(reader);
            if (acquired < 0) {
                return -3;
            }
            if (acquired == 0) {
                continue;
            }
            const char *input = (const char *)reader->span.data + reader->position;
            size_t available = reader->span.size - reader->position;
            const char *newline = memchr(input, '\n', available);
            size_t copied = newline ? (size_t)(newline - input) + 1 : available;
            if (append_carry(reader, input, copied) != 0) {
                return -3;
            }
            reader->position += copied;
            continue;
        }

        int acquired = acquire_span(reader);
        if (acquired < 0) {
            return -3;
        }
        if (acquired == 0) {
            return 0;
        }
        const char *input = (const char *)reader->span.data + reader->position;
        size_t available = reader->span.size - reader->position;
        const char *header = find_header(input, available);
        if (!header) {
            reader->position = reader->span.size;
            continue;
        }
        reader->position += (size_t)(header - input);
        input = header;
        available = reader->span.size - reader->position;
        record_layout layout;
        parse_status parsed = parse_layout(input, available, 0, &layout);
        if (parsed == PARSE_COMPLETE) {
            if (build_record(reader, bam, record, input, &layout) != 0) {
                return -3;
            }
            reader->position += layout.consumed;
            return 1;
        }
        if (parsed == PARSE_QUAL_ERROR) {
            reader->position += layout.consumed;
            return -2;
        }
        if (parsed == PARSE_QUAL_MISSING) {
            reader->position += layout.consumed;
            return -4;
        }
        if (parsed == PARSE_INVALID) {
            return -3;
        }
        if (append_carry(reader, input, available) != 0) {
            return -3;
        }
        reader->position = reader->span.size;
    }
}
