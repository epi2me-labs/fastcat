#include <string.h>
#include <stdlib.h>

#include <zlib.h>

#include "unity.h"
#include "htslib/kseq.h"

#ifndef KSEQ_DECLARED
#define KSEQ_DECLARED
KSEQ_DECLARE(gzFile)
#endif

#include "fc_read.h"

static fc_read* read_from_comment(const char* comment) {
    kseq_t seq = {0};
    seq.name.s = strdup("read1");
    seq.name.l = strlen(seq.name.s);
    seq.seq.s = strdup("ACGT");
    seq.seq.l = strlen(seq.seq.s);
    seq.qual.s = strdup("####"); // ASCII 35 -> phred 2 after offset
    seq.qual.l = strlen(seq.qual.s);
    seq.comment.s = strdup(comment);
    seq.comment.l = strlen(seq.comment.s);

    fc_read* fr = fc_read_from_kseq(&seq);

    free(seq.name.s);
    free(seq.seq.s);
    free(seq.qual.s);
    free(seq.comment.s);
    return fr;
}

static void test_kseq_to_bam_unmapped_aux(void) {
    fc_read* fr = read_from_comment("XY:Z:foo");
    bam1_t* bam = fr->bam;
    TEST_ASSERT_TRUE(bam->core.flag & BAM_FUNMAP);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.tid);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.pos);
    TEST_ASSERT_EQUAL_INT(4, bam->core.l_qseq);
    TEST_ASSERT_EQUAL_STRING("read1", bam_get_qname(bam));

    uint8_t* qual = bam_get_qual(bam);
    TEST_ASSERT_EQUAL_UINT8(2, qual[0]); // 35 - 33

    uint8_t* aux = bam_aux_get(bam, "XY");
    TEST_ASSERT_NOT_NULL(aux);
    const char* aux_str = bam_aux2Z(aux);
    TEST_ASSERT_NOT_NULL(aux_str);
    TEST_ASSERT_EQUAL_STRING("foo", aux_str);

    fc_read_destroy(fr);
}

static void test_fc_read_from_kseq_sam_comment_consistency(void) {
    fc_read* fr = read_from_comment("RG:Z:rg1\tSM:Z:barcode07\tpi:Z:parent-uuid\tML:B:C,1,2");

    uint8_t* aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("parent-uuid", bam_aux2Z(aux_pi));

    uint8_t* aux_ml = bam_aux_get(fr->bam, "ML");
    TEST_ASSERT_NOT_NULL(aux_ml);
    TEST_ASSERT_EQUAL_UINT32(2, bam_auxB_len(aux_ml));
    TEST_ASSERT_EQUAL_INT64(1, bam_auxB2i(aux_ml, 0));
    TEST_ASSERT_EQUAL_INT64(2, bam_auxB2i(aux_ml, 1));

    TEST_ASSERT_EQUAL_STRING("parent-uuid", fr->meta->parent_read);
    TEST_ASSERT_EQUAL_STRING("barcode07", fr->meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, fr->meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "pi:Z:parent-uuid"));
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "SM:Z:barcode07"));
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "ML:B:C,1,2"));

    fc_read_destroy(fr);
}

static void test_fc_read_from_kseq_key_value_fallback_consistency(void) {
    fc_read* fr = read_from_comment("runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model_version_id=rna004_130bps_fast@v5.1.0 parent_read_id=pr-001 barcode=barcode01 start_time=2023-01-01T00:00:00Z");

    uint8_t* aux_rg = bam_aux_get(fr->bam, "RG");
    TEST_ASSERT_NOT_NULL(aux_rg);
    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode01", bam_aux2Z(aux_rg));

    uint8_t* aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("pr-001", bam_aux2Z(aux_pi));

    uint8_t* aux_sm = bam_aux_get(fr->bam, "SM");
    TEST_ASSERT_NOT_NULL(aux_sm);
    TEST_ASSERT_EQUAL_STRING("barcode01", bam_aux2Z(aux_sm));

    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5", fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING("pr-001", fr->meta->parent_read);
    TEST_ASSERT_EQUAL_STRING("barcode01", fr->meta->barcode);
    TEST_ASSERT_EQUAL_INT(1, fr->meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", fr->meta->start_time);
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "pi:Z:pr-001"));
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "SM:Z:barcode01"));

    fc_read_destroy(fr);
}

static void test_fc_read_preserves_original_comment(void) {
    fc_read* fr = read_from_comment("raw\tcomment  with  spaces");
    TEST_ASSERT_NOT_NULL(fr->meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw\tcomment  with  spaces", fr->meta->orig_comment);

    fc_read_destroy(fr);
}

void bam_utils_suite(void) {
    RUN_TEST(test_kseq_to_bam_unmapped_aux);
    RUN_TEST(test_fc_read_from_kseq_sam_comment_consistency);
    RUN_TEST(test_fc_read_from_kseq_key_value_fallback_consistency);
    RUN_TEST(test_fc_read_preserves_original_comment);
}
