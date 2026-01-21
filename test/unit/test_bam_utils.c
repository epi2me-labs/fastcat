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
    const char* runid = "12345678-1234-1234-1234-123456789012";
    fc_read* fr = read_from_comment(
        "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1_barcode07\t"
        "SM:Z:barcode07\tpi:Z:parent-uuid\tML:B:C,1,2\t"
        "st:Z:2023-01-01T00:00:00Z\trn:i:7\tch:i:42\t"
        "PU:Z:FLO-MIN114\tLB:Z:sample123\tal:Z:alias1\tDS:Z:gpu:gpu0");

    uint8_t* aux_rg = bam_aux_get(fr->bam, "RG");
    TEST_ASSERT_NOT_NULL(aux_rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012_sup@v1_barcode07", bam_aux2Z(aux_rg));

    uint8_t* aux_sm = bam_aux_get(fr->bam, "SM");
    TEST_ASSERT_NOT_NULL(aux_sm);
    TEST_ASSERT_EQUAL_STRING("barcode07", bam_aux2Z(aux_sm));

    uint8_t* aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("parent-uuid", bam_aux2Z(aux_pi));

    uint8_t* aux_ml = bam_aux_get(fr->bam, "ML");
    TEST_ASSERT_NOT_NULL(aux_ml);
    TEST_ASSERT_EQUAL_UINT32(2, bam_auxB_len(aux_ml));
    TEST_ASSERT_EQUAL_INT64(1, bam_auxB2i(aux_ml, 0));
    TEST_ASSERT_EQUAL_INT64(2, bam_auxB2i(aux_ml, 1));

    uint8_t* aux_st = bam_aux_get(fr->bam, "st");
    TEST_ASSERT_NOT_NULL(aux_st);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", bam_aux2Z(aux_st));

    uint8_t* aux_rn = bam_aux_get(fr->bam, "rn");
    TEST_ASSERT_NOT_NULL(aux_rn);
    TEST_ASSERT_EQUAL_INT(7, bam_aux2i(aux_rn));

    uint8_t* aux_ch = bam_aux_get(fr->bam, "ch");
    TEST_ASSERT_NOT_NULL(aux_ch);
    TEST_ASSERT_EQUAL_INT(42, bam_aux2i(aux_ch));

    uint8_t* aux_pu = bam_aux_get(fr->bam, "PU");
    TEST_ASSERT_NOT_NULL(aux_pu);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", bam_aux2Z(aux_pu));

    uint8_t* aux_lb = bam_aux_get(fr->bam, "LB");
    TEST_ASSERT_NOT_NULL(aux_lb);
    TEST_ASSERT_EQUAL_STRING("sample123", bam_aux2Z(aux_lb));

    uint8_t* aux_al = bam_aux_get(fr->bam, "al");
    TEST_ASSERT_NOT_NULL(aux_al);
    TEST_ASSERT_EQUAL_STRING("alias1", bam_aux2Z(aux_al));

    uint8_t* aux_ds = bam_aux_get(fr->bam, "DS");
    TEST_ASSERT_NOT_NULL(aux_ds);
    TEST_ASSERT_EQUAL_STRING("gpu:gpu0", bam_aux2Z(aux_ds));

    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("parent-uuid", fr->meta->parent_read);
    TEST_ASSERT_EQUAL_STRING("barcode07", fr->meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, fr->meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", fr->meta->flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("alias1", fr->meta->barcode_alias);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", fr->meta->start_time);
    TEST_ASSERT_EQUAL_UINT64(7, fr->meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, fr->meta->channel);
    TEST_ASSERT_NOT_NULL(fr->meta->rg_info);
    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->rg_info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", fr->meta->rg_info->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode07", fr->meta->rg_info->barcode);
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "pi:Z:parent-uuid"));
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "SM:Z:barcode07"));
    TEST_ASSERT_NOT_NULL(strstr(fr->meta->tags_str->s, "ML:B:C,1,2"));

    fc_read_destroy(fr);
}

static void test_fc_read_from_kseq_key_value_fallback_consistency(void) {
    const char* runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5";
    const char* model = "rna004_130bps_fast@v5.1.0";
    fc_read* fr = read_from_comment(
        "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 "
        "basecall_model_version_id=rna004_130bps_fast@v5.1.0 "
        "flow_cell_id=FLO-MIN114 sample_id=sample123 "
        "barcode=barcode01 barcode_alias=alias1 read=7 ch=42 "
        "start_time=2023-01-01T00:00:00Z parent_read_id=pr-001 "
        "basecall_gpu=gpu0 empty= unknown=foo");

    uint8_t* aux_rg = bam_aux_get(fr->bam, "RG");
    TEST_ASSERT_NOT_NULL(aux_rg);
    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode01",
        bam_aux2Z(aux_rg));

    uint8_t* aux_pu = bam_aux_get(fr->bam, "PU");
    TEST_ASSERT_NOT_NULL(aux_pu);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", bam_aux2Z(aux_pu));

    uint8_t* aux_lb = bam_aux_get(fr->bam, "LB");
    TEST_ASSERT_NOT_NULL(aux_lb);
    TEST_ASSERT_EQUAL_STRING("sample123", bam_aux2Z(aux_lb));

    uint8_t* aux_al = bam_aux_get(fr->bam, "al");
    TEST_ASSERT_NOT_NULL(aux_al);
    TEST_ASSERT_EQUAL_STRING("alias1", bam_aux2Z(aux_al));

    uint8_t* aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("pr-001", bam_aux2Z(aux_pi));

    uint8_t* aux_sm = bam_aux_get(fr->bam, "SM");
    TEST_ASSERT_NOT_NULL(aux_sm);
    TEST_ASSERT_EQUAL_STRING("barcode01", bam_aux2Z(aux_sm));

    uint8_t* aux_rn = bam_aux_get(fr->bam, "rn");
    TEST_ASSERT_NOT_NULL(aux_rn);
    TEST_ASSERT_EQUAL_INT(7, bam_aux2i(aux_rn));

    uint8_t* aux_ch = bam_aux_get(fr->bam, "ch");
    TEST_ASSERT_NOT_NULL(aux_ch);
    TEST_ASSERT_EQUAL_INT(42, bam_aux2i(aux_ch));

    uint8_t* aux_st = bam_aux_get(fr->bam, "st");
    TEST_ASSERT_NOT_NULL(aux_st);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", bam_aux2Z(aux_st));

    uint8_t* aux_ds = bam_aux_get(fr->bam, "DS");
    TEST_ASSERT_NOT_NULL(aux_ds);
    TEST_ASSERT_EQUAL_STRING("gpu:gpu0", bam_aux2Z(aux_ds));

    uint8_t* aux_co = bam_aux_get(fr->bam, "CO");
    TEST_ASSERT_NOT_NULL(aux_co);
    TEST_ASSERT_EQUAL_STRING("empty= unknown=foo", bam_aux2Z(aux_co));

    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING(model, fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", fr->meta->flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("pr-001", fr->meta->parent_read);
    TEST_ASSERT_EQUAL_STRING("barcode01", fr->meta->barcode);
    TEST_ASSERT_EQUAL_INT(1, fr->meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", fr->meta->start_time);
    TEST_ASSERT_EQUAL_UINT64(7, fr->meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, fr->meta->channel);
    TEST_ASSERT_EQUAL_STRING("alias1", fr->meta->barcode_alias);
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
