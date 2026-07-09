#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "main.h"
#include "core/fc_fastqcomments.h"
#include "core/fc_read.h"
#include "core/fc_bam_utils.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

static void read_from_comment(fc_read *read, const char* comment) {
    if (!read->bam) {
        read->bam = bam_init1();
        TEST_ASSERT_NOT_NULL(read->bam);
    }
    const char *name = "dummy";
    const char *seq = "A";
    const char *qual_ascii = "!";
    TEST_ASSERT_TRUE(bam_set1(
        read->bam,
        strlen(name) + 1, name,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        strlen(seq), seq, qual_ascii,
        0) >= 0);

    uint8_t *qual = bam_get_qual(read->bam);
    qual[0] = (uint8_t)(qual[0] - 33);

    kstring_t c = { strlen(comment), strlen(comment), (char *)comment };
    if (c.l > 0) {
        TEST_ASSERT_EQUAL_INT(0, fastq_comment_to_bam(&c, read->bam));
    }

    fc_bam_ensure_mean_quality(read->bam, false);
    fc_read_reset_from_bam(read, NULL, NULL);

    char *tmp = realloc(read->meta->orig_comment, c.l + 1);
    TEST_ASSERT_NOT_NULL(tmp);
    read->meta->orig_comment = tmp;
    if (c.l > 0) {
        memcpy(read->meta->orig_comment, c.s, c.l);
    }
    read->meta->orig_comment[c.l] = '\0';
}

static void destroy_read(fc_read *read) {
    fc_read_destroy(read);
}

static const char *get_tags(fc_read *read) {
    const kstring_t *tags = fc_read_get_tags_str(read);
    TEST_ASSERT_NOT_NULL(tags);
    TEST_ASSERT_NOT_NULL(tags->s);
    return tags->s;
}

FC_TEST(fastqcomments_suite, test_parse_simple_comment)
{
    const char* runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5";
    const char* model = "rna004_130bps_fast@v5.1.0";
    char comment[256];
    FC_SNPRINTF_CHECK(comment, sizeof(comment),
        "runid=%s basecall_model_version_id=%s barcode=barcode12 start_time=2023-01-01T00:00:00Z du=1.25",
        runid, model);
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, comment);
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING(runid, meta->runid);
    TEST_ASSERT_EQUAL_STRING("barcode12", meta->barcode);
    TEST_ASSERT_EQUAL_INT(12, meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", meta->start_time);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.25, meta->duration_s);
    const char *tags = get_tags(read);
    TEST_ASSERT_NOT_NULL(strstr(tags, "RG:Z:e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode12"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "SM:Z:barcode12"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "st:Z:2023-01-01T00:00:00Z"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "du:f:1.25"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_original_comment_preserved)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, "raw_comment\twith\tspecials");
    read_meta *meta = read->meta;
    TEST_ASSERT_NOT_NULL(meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw_comment\twith\tspecials", meta->orig_comment);
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_unknown_key_goes_to_co)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model_version_id=rna004_130bps_fast@v5.1.0 unknown=foo");
    const char *tags = get_tags(read);
    TEST_ASSERT_NOT_NULL(strstr(tags, "RG:Z:e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "CO:Z:unknown=foo"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_empty_value_goes_to_co_and_doesnt_cause_hell_to_open)
{
    const char* runid = "12345678-1234-1234-1234-123456789012";
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read,
        "runid=12345678-1234-1234-1234-123456789012 empty= no_equals barcode=barcode07 "
        "fubar= read=7 unknown=foo ch=42 mux=3 du=2.5");
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING(runid, meta->runid);
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_EQUAL_UINT64(7, meta->read_number);
    TEST_ASSERT_EQUAL_UINT64(42, meta->channel);
    TEST_ASSERT_EQUAL_UINT64(3, meta->mux);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 2.5, meta->duration_s);

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

    uint8_t* mx_aux = bam_aux_get(read->bam, "mx");
    TEST_ASSERT_NOT_NULL(mx_aux);
    TEST_ASSERT_EQUAL_INT(3, bam_aux2i(mx_aux));

    uint8_t* du_aux = bam_aux_get(read->bam, "du");
    TEST_ASSERT_NOT_NULL(du_aux);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 2.5, bam_aux2f(du_aux));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_sam_style_tabs)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, "RG:Z:rg1\tSM:Z:barcode01\tZZ:Z:extra");
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("rg1", meta->rg);
    TEST_ASSERT_EQUAL_STRING("barcode01", meta->barcode);
    const char *tags = get_tags(read);
    TEST_ASSERT_NOT_NULL(strstr(tags, "RG:Z:rg1"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "SM:Z:barcode01"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "ZZ:Z:extra"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_barcode_unclassified)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, "barcode=unclassified");
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("unclassified", meta->barcode);
    TEST_ASSERT_EQUAL_INT(0, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(get_tags(read), "SM:Z:unclassified"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_barcode_empty_value_fallback)
{
    const char *runid = "12345678-1234-1234-1234-123456789012";
    const char *model = "sup@v1";
    char comment[256];
    FC_SNPRINTF_CHECK(comment, sizeof(comment),
        "runid=%s basecall_model_version_id=%s barcode=", runid, model);

    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, comment);
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("", meta->barcode);
    TEST_ASSERT_EQUAL_INT(0, meta->ibarcode);
    const char *tags = get_tags(read);
    TEST_ASSERT_NOT_NULL(strstr(tags, "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1"));
    TEST_ASSERT_NOT_NULL(strstr(tags, "SM:Z:"));

    uint8_t *sm_aux = bam_aux_get(read->bam, "SM");
    TEST_ASSERT_NOT_NULL(sm_aux);
    TEST_ASSERT_EQUAL_STRING("", bam_aux2Z(sm_aux));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_rg_populates_runid_and_basecaller)
{
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v1";
    char buf[256];
    FC_SNPRINTF_CHECK(buf, sizeof(buf), "RG:Z:%s", rg);
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, buf);
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", meta->basecaller);
    TEST_ASSERT_NULL(meta->rg_info->barcode);
    TEST_ASSERT_NOT_NULL(strstr(get_tags(read), "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_rg_fallback)
{
    const char* rg = "12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode07";
    char buf[256];
    FC_SNPRINTF_CHECK(buf, sizeof(buf), "RG:Z:%s\tSM:Z:barcode07", rg);
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, buf);
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012", meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v4.0.0", meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_NOT_NULL(strstr(get_tags(read), "RG:Z:12345678-1234-1234-1234-123456789012_sup@v4.0.0_barcode07"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_mixed_sam_and_key_value_prefers_sam)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read,
        "RG:Z:rg1\tSM:Z:barcode07\tparent_read_id=pr-001 start_time=2023-01-01T00:00:00Z");
    read_meta *meta = read->meta;
    TEST_ASSERT_EQUAL_STRING("barcode07", meta->barcode);
    TEST_ASSERT_EQUAL_INT(7, meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("", meta->start_time);
    TEST_ASSERT_NULL(bam_aux_get(read->bam, "pi"));
    TEST_ASSERT_NULL(strstr(get_tags(read), "CO:Z:"));
    destroy_read(read);
}

FC_TEST(fastqcomments_suite, test_reheader_roundtrip_preserves_aux)
{
    const char* comment =
        "runid=run-123 basecall_model_version_id=model@v1 barcode=barcode07 "
        "start_time=2023-01-01T00:00:00Z parent_read_id=pr-001 read=7 ch=42 "
        "flow_cell_id=FLO-MIN114 basecall_gpu=gpu0 barcode_alias=alias1 unknown=foo";
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);
    read_from_comment(read, comment);
    const kstring_t *tags = fc_read_get_tags_str(read);
    TEST_ASSERT_NOT_NULL(tags);
    TEST_ASSERT_TRUE(tags->l > 0);

    kstring_t reheader = {0, 0, NULL};
    ksprintf(&reheader, "%s", tags->s);

    bam1_t* b2 = bam_init1();
    TEST_ASSERT_NOT_NULL(b2);
    const char* qname = "dummy2";
    const char* seq = "A";
    const char* qual = "!";
    TEST_ASSERT_TRUE(bam_set1(b2, strlen(qname) + 1, qname,
        BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, strlen(seq), seq, qual, 0) >= 0);
    TEST_ASSERT_EQUAL_INT(0, fastq_comment_to_bam(&reheader, b2));

    // compare the lazily materialized SAM-style aux tag strings used by fastcat
    fc_read *read2 = fc_read_init();
    TEST_ASSERT_NOT_NULL(read2);
    bam_destroy1(read2->bam);
    read2->bam = b2;
    fc_bam_ensure_mean_quality(read2->bam, false);
    fc_read_reset_from_bam(read2, NULL, NULL);
    const kstring_t *aux1 = fc_read_get_tags_str(read);
    const kstring_t *aux2 = fc_read_get_tags_str(read2);
    TEST_ASSERT_NOT_NULL(aux1);
    TEST_ASSERT_NOT_NULL(aux2);
    TEST_ASSERT_EQUAL_STRING(aux1->s, aux2->s);

    ks_free(&reheader);
    read2->bam = NULL;
    bam_destroy1(b2);
    fc_read_destroy(read2);
    destroy_read(read);
}
