#include <string.h>
#include <stdlib.h>

#include <zlib.h>

#include "unity.h"
#include "main.h"
#include "htslib/kseq.h"

#ifndef KSEQ_DECLARED
#define KSEQ_DECLARED
KSEQ_DECLARE(gzFile)
#endif

#include "core/fc_read.h"
#include "core/fc_fastqcomments.h"
#include "core/fc_bam_utils.h"

static void read_from_comment(fc_read *fr, const char* comment) {
    if (!fr->bam) {
        fr->bam = bam_init1();
        TEST_ASSERT_NOT_NULL(fr->bam);
    }
    const char *name = "read1";
    const char *seq = "ACGT";
    const char *qual_ascii = "####"; // ASCII 35 -> phred 2 after offset
    TEST_ASSERT_TRUE(bam_set1(
        fr->bam,
        strlen(name) + 1, name,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        strlen(seq), seq, qual_ascii,
        0) >= 0);

    uint8_t *qual = bam_get_qual(fr->bam);
    for (int i = 0; i < fr->bam->core.l_qseq; ++i) {
        qual[i] = (uint8_t)(qual[i] - 33);
    }

    kstring_t c = { strlen(comment), strlen(comment), (char *)comment };
    if (c.l > 0) {
        TEST_ASSERT_EQUAL_INT(0, fastq_comment_to_bam(&c, fr->bam));
    }

    fc_bam_ensure_mean_quality(fr->bam, false);
    fc_read_reset_from_bam(fr, NULL, NULL);

    char *tmp = realloc(fr->meta->orig_comment, c.l + 1);
    TEST_ASSERT_NOT_NULL(tmp);
    fr->meta->orig_comment = tmp;
    if (c.l > 0) {
        memcpy(fr->meta->orig_comment, c.s, c.l);
    }
    fr->meta->orig_comment[c.l] = '\0';
}

static void destroy_read(fc_read *fr) {
    fc_read_destroy(fr);
}

static const char *get_tags(fc_read *read) {
    const kstring_t *tags = fc_read_get_tags_str(read);
    TEST_ASSERT_NOT_NULL(tags);
    TEST_ASSERT_NOT_NULL(tags->s);
    return tags->s;
}

FC_TEST(bam_utils_suite, test_kseq_to_bam_unmapped_aux) {
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "XY:Z:foo");
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

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_sam_comment_consistency) {
    const char* runid = "12345678-1234-1234-1234-123456789012";
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr,
        "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1_barcode07\t"
        "SM:Z:barcode07\tpi:Z:parent-uuid\tML:B:C,1,2\t"
        "st:Z:2023-01-01T00:00:00Z\tdu:f:1.25\trn:i:7\tch:i:42\t"
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

    uint8_t* aux_du = bam_aux_get(fr->bam, "du");
    TEST_ASSERT_NOT_NULL(aux_du);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, bam_aux2f(aux_du));

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
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.25, fr->meta->duration_s);
    TEST_ASSERT_EQUAL_UINT64(7, fr->meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, fr->meta->channel);
    TEST_ASSERT_NOT_NULL(fr->meta->rg_info);
    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->rg_info->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", fr->meta->rg_info->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode07", fr->meta->rg_info->barcode);
    const char *tags = get_tags(fr);
    TEST_ASSERT_NOT_NULL(strstr(tags, "pi:Z:parent-uuid"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "SM:Z:barcode07"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "ML:B:C,1,2"));

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_key_value_fallback_consistency) {
    const char* runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5";
    const char* model = "rna004_130bps_fast@v5.1.0";
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr,
        "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 "
        "basecall_model_version_id=rna004_130bps_fast@v5.1.0 "
        "flow_cell_id=FLO-MIN114 sample_id=sample123 "
        "barcode=barcode01 barcode_alias=alias1 read=7 ch=42 "
        "start_time=2023-01-01T00:00:00Z du=2.5 parent_read_id=pr-001 "
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

    uint8_t* aux_du = bam_aux_get(fr->bam, "du");
    TEST_ASSERT_NOT_NULL(aux_du);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 2.5, bam_aux2f(aux_du));

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
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 2.5, fr->meta->duration_s);
    TEST_ASSERT_EQUAL_UINT64(7, fr->meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, fr->meta->channel);
    TEST_ASSERT_EQUAL_STRING("alias1", fr->meta->barcode_alias);
    const char *tags = get_tags(fr);
    TEST_ASSERT_NOT_NULL(strstr(tags, "pi:Z:pr-001"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "SM:Z:barcode01"));

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fc_read_preserves_original_comment) {
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "raw\tcomment  with  spaces");
    TEST_ASSERT_NOT_NULL(fr->meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw\tcomment  with  spaces", fr->meta->orig_comment);

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_reuses_preallocated_bam) {
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "XY:Z:foo");
    bam1_t *first_bam = fr->bam;
    TEST_ASSERT_NOT_NULL(first_bam);

    read_from_comment(fr, "XY:Z:bar");
    TEST_ASSERT_EQUAL_PTR(first_bam, fr->bam);

    uint8_t* aux = bam_aux_get(fr->bam, "XY");
    TEST_ASSERT_NOT_NULL(aux);
    TEST_ASSERT_EQUAL_STRING("bar", bam_aux2Z(aux));
    destroy_read(fr);
}
