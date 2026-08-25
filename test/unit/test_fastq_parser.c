#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "core/fc_fastq_parser.h"
#include "core/kseq_gz.h"
#include "main.h"
#include "test_fs_utils.h"
#include "unity.h"

static void write_fastq_gzip_bytes(
    const char *contents, size_t length, char *path, size_t path_capacity)
{
    char temporary[FC_TEST_TEMP_PATH_CAPACITY] = "/tmp/fastcat_parserXXXXXX";
    int fd = fc_test_mkstemp(temporary, sizeof(temporary));
    TEST_ASSERT_NOT_EQUAL(-1, fd);
    close(fd);

    gzFile gzip = gzopen(temporary, "wb");
    TEST_ASSERT_NOT_NULL(gzip);
    TEST_ASSERT_EQUAL_INT((int)length, gzwrite(gzip, contents, (unsigned)length));
    TEST_ASSERT_EQUAL_INT(Z_OK, gzclose(gzip));
    FC_SNPRINTF_CHECK(path, path_capacity, "%s.fastq.gz", temporary);
    TEST_ASSERT_EQUAL_INT(0, rename(temporary, path));
}

static void write_fastq_gzip(const char *contents, char *path, size_t path_capacity)
{
    write_fastq_gzip_bytes(contents, strlen(contents), path, path_capacity);
}

static void assert_bam_sequence(const bam1_t *bam, const char *expected)
{
    const uint8_t *packed = bam_get_seq(bam);
    size_t length = strlen(expected);
    TEST_ASSERT_EQUAL_UINT64(length, bam->core.l_qseq);
    for (size_t i = 0; i < length; ++i) {
        TEST_ASSERT_EQUAL_CHAR(expected[i], seq_nt16_str[bam_seqi(packed, i)]);
    }
}

typedef enum {
    FASTQ_EVENT_RECORD,
    FASTQ_EVENT_EOF,
    FASTQ_EVENT_QUALITY_MISSING,
    FASTQ_EVENT_QUALITY_TRUNCATED,
    FASTQ_EVENT_MALFORMED,
    FASTQ_EVENT_OVERFLOW
} fastq_event;

static fastq_event next_kseq_event(kseq_t *reader, bam1_t *bam)
{
    int status = kseq_read(reader);
    if (status == -1) {
        return FASTQ_EVENT_EOF;
    }
    if (status == -2) {
        return FASTQ_EVENT_QUALITY_TRUNCATED;
    }
    if (status == -3) {
        return FASTQ_EVENT_MALFORMED;
    }
    if (status == -4) {
        return FASTQ_EVENT_OVERFLOW;
    }
    if (reader->qual.l == 0) {
        return FASTQ_EVENT_QUALITY_MISSING;
    }
    return fc_fastq_record_to_bam(bam,
               reader->name.s, reader->name.l,
               reader->comment.s, reader->comment.l,
               reader->seq.s, reader->qual.s, reader->seq.l)
            < 0
        ? FASTQ_EVENT_MALFORMED
        : FASTQ_EVENT_RECORD;
}

static fastq_event next_parser_event(fc_fastq_parser *reader,
    bam1_t *bam, fc_fastq_record *record)
{
    int status = fc_fastq_parser_next(reader, bam, record);
    switch (status) {
    case 1:
        return FASTQ_EVENT_RECORD;
    case 0:
        return FASTQ_EVENT_EOF;
    case -2:
        return FASTQ_EVENT_QUALITY_TRUNCATED;
    case -4:
        return FASTQ_EVENT_QUALITY_MISSING;
    default:
        return FASTQ_EVENT_MALFORMED;
    }
}

static void assert_parser_matches_kseq_path(const char *path)
{
    fc_fastq_stream *kseq_stream = fc_fastq_stream_open_path(path, 2);
    fc_fastq_stream *parser_stream = fc_fastq_stream_open_path(path, 2);
    TEST_ASSERT_NOT_NULL(kseq_stream);
    TEST_ASSERT_NOT_NULL(parser_stream);
    kseq_t *kseq = kseq_init(kseq_stream);
    fc_fastq_parser *parser = fc_fastq_parser_create(parser_stream);
    TEST_ASSERT_NOT_NULL(kseq);
    TEST_ASSERT_NOT_NULL(parser);
    bam1_t *kseq_bam = bam_init1();
    bam1_t *parser_bam = bam_init1();
    TEST_ASSERT_NOT_NULL(kseq_bam);
    TEST_ASSERT_NOT_NULL(parser_bam);

    for (size_t event_index = 0; event_index < 32; ++event_index) {
        fc_fastq_record parsed_record = { 0 };
        fastq_event expected = next_kseq_event(kseq, kseq_bam);
        fastq_event actual = next_parser_event(parser, parser_bam, &parsed_record);
        TEST_ASSERT_EQUAL_INT(expected, actual);
        if (expected == FASTQ_EVENT_RECORD) {
            TEST_ASSERT_EQUAL_UINT64(kseq->name.l, parsed_record.name_len);
            TEST_ASSERT_EQUAL_UINT64(kseq->comment.l, parsed_record.comment_len);
            TEST_ASSERT_EQUAL_UINT64(kseq->seq.l, parsed_record.sequence_len);
            if (kseq->name.l != 0) {
                TEST_ASSERT_EQUAL_MEMORY(kseq->name.s, parsed_record.name, kseq->name.l);
            }
            if (kseq->comment.l != 0) {
                TEST_ASSERT_EQUAL_MEMORY(kseq->comment.s, parsed_record.comment, kseq->comment.l);
            }
            if (kseq->seq.l != 0) {
                TEST_ASSERT_EQUAL_MEMORY(kseq->seq.s, parsed_record.sequence, kseq->seq.l);
                TEST_ASSERT_EQUAL_MEMORY(kseq->qual.s, parsed_record.quality, kseq->qual.l);
            }
            TEST_ASSERT_EQUAL_MEMORY(
                &kseq_bam->core, &parser_bam->core, sizeof(kseq_bam->core));
            TEST_ASSERT_EQUAL_INT(kseq_bam->l_data, parser_bam->l_data);
            TEST_ASSERT_EQUAL_MEMORY(kseq_bam->data, parser_bam->data, kseq_bam->l_data);
        }
        if (expected == FASTQ_EVENT_EOF || expected == FASTQ_EVENT_MALFORMED
            || expected == FASTQ_EVENT_OVERFLOW) {
            break;
        }
        TEST_ASSERT_LESS_THAN_size_t(31, event_index);
    }

    bam_destroy1(parser_bam);
    bam_destroy1(kseq_bam);
    fc_fastq_parser_destroy(parser);
    kseq_destroy(kseq);
    int direct_close_status = fc_fastq_stream_close(parser_stream);
    int kseq_close_status = fc_fastq_stream_close(kseq_stream);
    TEST_ASSERT_EQUAL_INT(kseq_close_status, direct_close_status);
}

static void assert_parser_matches_kseq(const char *contents)
{
    char path[FC_TEST_TEMP_PATH_CAPACITY];
    write_fastq_gzip(contents, path, sizeof(path));
    assert_parser_matches_kseq_path(path);
    TEST_ASSERT_EQUAL_INT(0, remove(path));
}

FC_TEST(fastq_parser_suite, test_fastq_parser_builds_bam_from_borrowed_spans)
{
    static const char input[] = "@read1 runid=abc ch=7\n"
                                "ACGTNRY\n"
                                "+\n"
                                "!\"#$%&'\n"
                                "@read2\tRG:Z:group\r\n"
                                "ACG\r\n"
                                "TN\r\n"
                                "+\r\n"
                                "III\r\n"
                                "II\r\n";
    char path[FC_TEST_TEMP_PATH_CAPACITY];
    write_fastq_gzip(input, path, sizeof(path));

    fc_fastq_stream *stream = fc_fastq_stream_open_path(path, 2);
    TEST_ASSERT_NOT_NULL(stream);
    fc_fastq_parser *parser = fc_fastq_parser_create(stream);
    TEST_ASSERT_NOT_NULL(parser);
    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    fc_fastq_record record = { 0 };

    TEST_ASSERT_EQUAL_INT(1, fc_fastq_parser_next(parser, bam, &record));
    TEST_ASSERT_EQUAL_UINT64(5, record.name_len);
    TEST_ASSERT_EQUAL_MEMORY("read1", record.name, record.name_len);
    TEST_ASSERT_EQUAL_MEMORY("ACGTNRY", record.sequence, record.sequence_len);
    TEST_ASSERT_EQUAL_MEMORY("!\"#$%&'", record.quality, record.sequence_len);
    TEST_ASSERT_EQUAL_STRING("read1", bam_get_qname(bam));
    assert_bam_sequence(bam, "ACGTNRY");
    TEST_ASSERT_EQUAL_UINT8(0, bam_get_qual(bam)[0]);
    TEST_ASSERT_EQUAL_UINT8(6, bam_get_qual(bam)[6]);

    TEST_ASSERT_EQUAL_INT(1, fc_fastq_parser_next(parser, bam, &record));
    TEST_ASSERT_EQUAL_MEMORY("read2", record.name, record.name_len);
    TEST_ASSERT_EQUAL_UINT64(5, record.sequence_len);
    TEST_ASSERT_EQUAL_MEMORY("ACGTN", record.sequence, record.sequence_len);
    TEST_ASSERT_EQUAL_MEMORY("IIIII", record.quality, record.sequence_len);
    TEST_ASSERT_EQUAL_STRING("read2", bam_get_qname(bam));
    assert_bam_sequence(bam, "ACGTN");
    TEST_ASSERT_NOT_NULL(bam_aux_get(bam, "RG"));

    TEST_ASSERT_EQUAL_INT(0, fc_fastq_parser_next(parser, bam, &record));
    bam_destroy1(bam);
    fc_fastq_parser_destroy(parser);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_stream_close(stream));
    TEST_ASSERT_EQUAL_INT(0, remove(path));
}

FC_TEST(fastq_parser_suite, test_fastq_parser_carries_record_across_views)
{
    const size_t sequence_length = 600000;
    const size_t input_length = 6 + sequence_length + 3 + sequence_length + 1;
    char *input = malloc(input_length);
    TEST_ASSERT_NOT_NULL(input);
    size_t offset = 0;
    memcpy(input + offset, "@long\n", 6);
    offset += 6;
    memset(input + offset, 'A', sequence_length);
    offset += sequence_length;
    memcpy(input + offset, "\n+\n", 3);
    offset += 3;
    memset(input + offset, 'I', sequence_length);
    offset += sequence_length;
    input[offset++] = '\n';
    TEST_ASSERT_EQUAL_UINT64(input_length, offset);

    char path[FC_TEST_TEMP_PATH_CAPACITY];
    write_fastq_gzip_bytes(input, input_length, path, sizeof(path));
    free(input);

    // One worker selects Purrflate's direct zlib-ng path, whose view is 1 MiB.
    fc_fastq_stream *stream = fc_fastq_stream_open_path(path, 1);
    TEST_ASSERT_NOT_NULL(stream);
    fc_fastq_parser *parser = fc_fastq_parser_create(stream);
    TEST_ASSERT_NOT_NULL(parser);
    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    fc_fastq_record record = { 0 };

    TEST_ASSERT_EQUAL_INT(1, fc_fastq_parser_next(parser, bam, &record));
    TEST_ASSERT_EQUAL_UINT64(sequence_length, record.sequence_len);
    TEST_ASSERT_EQUAL_CHAR('A', record.sequence[0]);
    TEST_ASSERT_EQUAL_CHAR('A', record.sequence[sequence_length - 1]);
    TEST_ASSERT_EQUAL_CHAR('I', record.quality[0]);
    TEST_ASSERT_EQUAL_CHAR('I', record.quality[sequence_length - 1]);
    TEST_ASSERT_EQUAL_UINT8(1, bam_seqi(bam_get_seq(bam), 0));
    TEST_ASSERT_EQUAL_UINT8(40, bam_get_qual(bam)[sequence_length - 1]);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_parser_next(parser, bam, &record));

    bam_destroy1(bam);
    fc_fastq_parser_destroy(parser);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_stream_close(stream));
    TEST_ASSERT_EQUAL_INT(0, remove(path));
}

FC_TEST(fastq_parser_suite, test_fastq_parser_matches_kseq_multiline_and_crlf)
{
    assert_parser_matches_kseq(
        "@multi comment with spaces\r\n"
        "AC\r\n"
        "\r\n"
        "GT\r\n"
        "+different identifier is permitted\r\n"
        "!@\r\n"
        "+I\r\n"
        "@final\n"
        "N\n"
        "+\n"
        "#");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_matches_kseq_header_like_quality)
{
    assert_parser_matches_kseq(
        "@symbols\n"
        "ACGT\n"
        "+ignored plus comment\n"
        "@+ I\n"
        "@next\n"
        "TGCA\n"
        "+\n"
        "IIII\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_matches_kseq_empty_read)
{
    assert_parser_matches_kseq(
        "@empty\n"
        "\n"
        "+\n"
        "\n"
        "@valid\n"
        "A\n"
        "+\n"
        "I\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_matches_kseq_empty_input)
{
    assert_parser_matches_kseq("");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_skips_fasta_record_like_kseq)
{
    assert_parser_matches_kseq(
        ">fasta comment\n"
        "AC\n"
        "GT\n"
        "@fastq\n"
        "A\n"
        "+\n"
        "I\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_matches_kseq_empty_name)
{
    assert_parser_matches_kseq(
        "@ comment only\n"
        "A\n"
        "+\n"
        "I\n"
        "@valid\n"
        "C\n"
        "+\n"
        "I\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_reports_header_only_eof_like_kseq)
{
    assert_parser_matches_kseq("@header-only\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_recovers_from_leading_noise_like_kseq)
{
    assert_parser_matches_kseq(
        "blank and damaged preamble\n"
        "\n"
        "@valid\n"
        "ACGT\n"
        "+\n"
        "IIII\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_recovers_missing_separator_like_kseq)
{
    assert_parser_matches_kseq(
        "@missing-plus\n"
        "ACGT\n"
        "@recovered\n"
        "TGCA\n"
        "+\n"
        "IIII\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_reports_missing_separator_at_eof_like_kseq)
{
    assert_parser_matches_kseq(
        "@missing-plus\n"
        "ACGT\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_reports_absent_quality_at_eof_like_kseq)
{
    assert_parser_matches_kseq(
        "@absent-quality\n"
        "ACGT\n"
        "+\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_reports_short_quality_at_eof_like_kseq)
{
    assert_parser_matches_kseq(
        "@short-quality\n"
        "ACGT\n"
        "+\n"
        "II\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_recovers_after_long_quality_like_kseq)
{
    assert_parser_matches_kseq(
        "@long-quality\n"
        "ACGT\n"
        "+\n"
        "IIIII\n"
        "@recovered\n"
        "TGCA\n"
        "+\n"
        "IIII\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_resynchronises_after_short_quality_like_kseq)
{
    assert_parser_matches_kseq(
        "@short-quality\n"
        "ACGT\n"
        "+\n"
        "II\n"
        "@consumed-as-quality\n"
        "TGCA\n"
        "+\n"
        "IIII\n"
        "@recovered\n"
        "N\n"
        "+\n"
        "I\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_recovers_from_inter_record_noise_like_kseq)
{
    assert_parser_matches_kseq(
        "@first\n"
        "A\n"
        "+\n"
        "I\n"
        "not a header\n"
        "another damaged line\n"
        "@recovered\n"
        "C\n"
        "+\n"
        "I\n");
}

FC_TEST(fastq_parser_suite, test_fastq_parser_reports_corrupt_gzip_like_kseq)
{
    static const char input[] = "@read\nACGT\n+\nIIII\n";
    char path[FC_TEST_TEMP_PATH_CAPACITY];
    struct stat file_info;
    write_fastq_gzip(input, path, sizeof(path));
    TEST_ASSERT_EQUAL_INT(0, stat(path, &file_info));
    TEST_ASSERT_GREATER_THAN_INT(8, file_info.st_size);
    TEST_ASSERT_EQUAL_INT(0, truncate(path, file_info.st_size - 8));
    assert_parser_matches_kseq_path(path);
    TEST_ASSERT_EQUAL_INT(0, remove(path));
}
