#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "fastqcomments.h"
#include "fc_read.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

static fc_read* read_from_comment(const char* comment) {
    kseq_t seq = {0};
    seq.name.s = strdup("dummy");
    seq.name.l = strlen(seq.name.s);
    seq.seq.s = strdup("A");
    seq.seq.l = strlen(seq.seq.s);
    seq.qual.s = strdup("!");
    seq.qual.l = strlen(seq.qual.s);
    seq.comment.s = strdup(comment);
    seq.comment.l = strlen(seq.comment.s);

    fc_read* read = fc_read_from_kseq(&seq);

    free(seq.name.s);
    free(seq.seq.s);
    free(seq.qual.s);
    free(seq.comment.s);
    return read;
}

static void test_parse_simple_comment(void) {
    const char* runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5";
    const char* model = "rna004_130bps_fast@v5.1.0";
    char comment[256];
    snprintf(comment, sizeof(comment),
        "runid=%s basecall_model_version_id=%s barcode=barcode12 start_time=2023-01-01T00:00:00Z",
        runid, model);
    fc_read* read = read_from_comment(comment);
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING(runid, meta->runid);
    TEST_ASSERT_EQUAL_STRING("barcode12", meta->barcode);
    TEST_ASSERT_EQUAL_INT(12, meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", meta->start_time);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode12"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "SM:Z:barcode12"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "st:Z:2023-01-01T00:00:00Z"));
    fc_read_destroy(read);
}

static void test_original_comment_preserved(void) {
    fc_read* read = read_from_comment("raw_comment\twith\tspecials");
    read_meta meta = read->meta;
    TEST_ASSERT_NOT_NULL(meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw_comment\twith\tspecials", meta->orig_comment);
    fc_read_destroy(read);
}

static void test_unknown_key_goes_to_co(void) {
    fc_read* read = read_from_comment("runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model_version_id=rna004_130bps_fast@v5.1.0 unknown=foo");
    read_meta meta = read->meta;
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "CO:Z:unknown=foo"));
    fc_read_destroy(read);
}

static void test_empty_value_goes_to_co_and_doesnt_cause_hell_to_open(void) {
    const char* runid = "12345678-1234-1234-1234-123456789012";
    fc_read* read = read_from_comment(
        "runid=12345678-1234-1234-1234-123456789012 empty= no_equals barcode=barcode07 "
        "fubar= read=7 unknown=foo ch=42");
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING(runid, meta->runid);
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_EQUAL_UINT64(7, meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, meta->channel);

    uint8_t* co_aux = bam_aux_get(read->bam, "CO");
    TEST_ASSERT_NOT_NULL(co_aux);
    TEST_ASSERT_EQUAL_STRING("empty= no_equals fubar= unknown=foo", bam_aux2Z(co_aux));

    uint8_t* sm_aux = bam_aux_get(read->bam, "SM");
    TEST_ASSERT_NOT_NULL(sm_aux);
    TEST_ASSERT_EQUAL_STRING("barcode07", bam_aux2Z(sm_aux));

    uint8_t* rn_aux = bam_aux_get(read->bam, "rn");
    TEST_ASSERT_NOT_NULL(rn_aux);
    TEST_ASSERT_EQUAL_INT(7, bam_aux2i(rn_aux));

    uint8_t* ch_aux = bam_aux_get(read->bam, "ch");
    TEST_ASSERT_NOT_NULL(ch_aux);
    TEST_ASSERT_EQUAL_INT(42, bam_aux2i(ch_aux));
    fc_read_destroy(read);
}

static void test_sam_style_tabs(void) {
    fc_read* read = read_from_comment("RG:Z:rg1\tSM:Z:barcode01\tZZ:Z:extra");
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("rg1", meta->rg);
    TEST_ASSERT_EQUAL_STRING("barcode01", meta->barcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:rg1"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "SM:Z:barcode01"));
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "ZZ:Z:extra"));
    fc_read_destroy(read);
}

static void test_barcode_unclassified(void) {
    fc_read* read = read_from_comment("barcode=unclassified");
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("unclassified", meta->barcode);
    TEST_ASSERT_EQUAL_INT(0, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "SM:Z:unclassified"));
    fc_read_destroy(read);
}

static void test_rg_populates_runid_and_basecaller(void) {
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v1";
    char buf[256];
    snprintf(buf, sizeof(buf), "RG:Z:%s", rg);
    fc_read* read = read_from_comment(buf);
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", meta->basecaller);
    TEST_ASSERT_NULL(meta->rg_info->barcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1"));
    fc_read_destroy(read);
}

static void test_rg_fallback(void) {
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode07";
    char buf[256];
    snprintf(buf, sizeof(buf), "RG:Z:%s\tSM:Z:barcode07", rg);
    fc_read* read = read_from_comment(buf);
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0", meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(meta->tags_str->s, "RG:Z:12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode07"));
    fc_read_destroy(read);
}

static void test_mixed_sam_and_key_value_prefers_sam(void) {
    fc_read* read = read_from_comment(
        "RG:Z:rg1\tSM:Z:barcode07\tparent_read_id=pr-001 start_time=2023-01-01T00:00:00Z");
    read_meta meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("", meta->start_time);
    TEST_ASSERT_NULL(bam_aux_get(read->bam, "pi"));
    TEST_ASSERT_NULL(strstr(meta->tags_str->s, "CO:Z:"));
    fc_read_destroy(read);
}

static void test_reheader_roundtrip_preserves_aux(void) {
    const char* comment =
        "runid=run-123 basecall_model_version_id=model@v1 barcode=barcode07 "
        "start_time=2023-01-01T00:00:00Z parent_read_id=pr-001 read=7 ch=42 "
        "flow_cell_id=FLO-MIN114 basecall_gpu=gpu0 barcode_alias=alias1 unknown=foo";
    fc_read* read = read_from_comment(comment);
    TEST_ASSERT_NOT_NULL(read->meta->tags_str);
    TEST_ASSERT_TRUE(read->meta->tags_str->l > 0);

    kstring_t reheader = {0, 0, NULL};
    ksprintf(&reheader, "%s", read->meta->tags_str->s);

    bam1_t* b2 = bam_init1();
    TEST_ASSERT_NOT_NULL(b2);
    const char* qname = "dummy2";
    const char* seq = "A";
    const char* qual = "!";
    TEST_ASSERT_TRUE(bam_set1(b2, strlen(qname) + 1, qname,
        BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, strlen(seq), seq, qual, 0) >= 0);
    TEST_ASSERT_EQUAL_INT(0, fastq_comment_to_bam(&reheader, b2));

    kstring_t aux1 = {0, 0, NULL};
    kstring_t aux2 = {0, 0, NULL};
    TEST_ASSERT_TRUE(sam_format_aux(read->bam, &aux1) >= 0);
    TEST_ASSERT_TRUE(sam_format_aux(b2, &aux2) >= 0);
    TEST_ASSERT_EQUAL_STRING(aux1.s, aux2.s);

    ks_free(&aux1);
    ks_free(&aux2);
    ks_free(&reheader);
    bam_destroy1(b2);
    fc_read_destroy(read);
}

void fastqcomments_suite(void) {
    RUN_TEST(test_parse_simple_comment);
    RUN_TEST(test_original_comment_preserved);
    RUN_TEST(test_unknown_key_goes_to_co);
    RUN_TEST(test_empty_value_goes_to_co_and_doesnt_cause_hell_to_open);
    RUN_TEST(test_sam_style_tabs);
    RUN_TEST(test_barcode_unclassified);
    RUN_TEST(test_rg_populates_runid_and_basecaller);
    RUN_TEST(test_rg_fallback);
    RUN_TEST(test_mixed_sam_and_key_value_prefers_sam);
    RUN_TEST(test_reheader_roundtrip_preserves_aux);
}
