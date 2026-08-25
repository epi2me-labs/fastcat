#include <string.h>
#include <stdlib.h>

#include "unity.h"
#include "main.h"

#include "core/fc_bam_rg_cache.h"
#include "core/fc_bam_utils.h"
#include "core/fc_fastqcomments.h"
#include "core/fc_fastq_parser.h"
#include "core/fc_read.h"

static void read_from_comment(fc_read *fr, const char* comment)
{
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
    fc_read_reset_from_bam(fr, NULL, fc_test_rg_cache(), &c);

    char *tmp = realloc(fr->meta->orig_comment, c.l + 1);
    TEST_ASSERT_NOT_NULL(tmp);
    fr->meta->orig_comment = tmp;
    if (c.l > 0) {
        memcpy(fr->meta->orig_comment, c.s, c.l);
    }
    fr->meta->orig_comment[c.l] = '\0';
}

static void destroy_read(fc_read *fr)
{
    fc_read_destroy(fr);
}

static sam_hdr_t *header_from_text(const char *text)
{
    sam_hdr_t *header = sam_hdr_parse((int)strlen(text), text);
    TEST_ASSERT_NOT_NULL(header);
    return header;
}

static const char *get_tags(fc_read *read)
{
    const kstring_t *tags = fc_read_get_tags_str(read);
    TEST_ASSERT_NOT_NULL(tags);
    TEST_ASSERT_NOT_NULL(tags->s);
    return tags->s;
}

static fc_read *read_from_rg_header_ds(const char *rg_id, const char *ds_text, bam_hdr_t **hdr_out, fc_bam_rg_cache **cache_out)
{
    char hdr_text[1024];
    FC_SNPRINTF_CHECK(
        hdr_text,
        sizeof(hdr_text),
        "@HD\tVN:1.6\n@RG\tID:%s\tDS:%s\n",
        rg_id,
        ds_text);

    bam_hdr_t *hdr = sam_hdr_parse((int)strlen(hdr_text), hdr_text);
    TEST_ASSERT_NOT_NULL(hdr);

    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    const char *name = "read1";
    const char *seq = "ACGT";
    const char *qual_ascii = "####";
    TEST_ASSERT_TRUE(bam_set1(
        bam,
        strlen(name) + 1, name,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        strlen(seq), seq, qual_ascii,
        0) >= 0);
    uint8_t *qual = bam_get_qual(bam);
    for (int i = 0; i < bam->core.l_qseq; ++i) {
        qual[i] = (uint8_t)(qual[i] - 33);
    }
    TEST_ASSERT_TRUE(bam_aux_update_str(bam, "RG", (int)strlen(rg_id) + 1, rg_id) >= 0);

    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    bam_destroy1(fr->bam);
    fr->bam = bam;

    fc_bam_rg_cache *cache = fc_bam_rg_cache_create(hdr);
    TEST_ASSERT_NOT_NULL(cache);
    fc_bam_ensure_mean_quality(fr->bam, false);
    fc_read_reset_from_bam(fr, hdr, cache, NULL);
    *hdr_out = hdr;
    *cache_out = cache;
    return fr;
}

FC_TEST(bam_utils_suite, test_unmapped_bam_aux)
{
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "XY:Z:foo");
    bam1_t *bam = fr->bam;
    TEST_ASSERT_TRUE(bam->core.flag & BAM_FUNMAP);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.tid);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.pos);
    TEST_ASSERT_EQUAL_INT(4, bam->core.l_qseq);
    TEST_ASSERT_EQUAL_STRING("read1", bam_get_qname(bam));

    uint8_t *qual = bam_get_qual(bam);
    TEST_ASSERT_EQUAL_UINT8(2, qual[0]); // 35 - 33

    uint8_t *aux = bam_aux_get(bam, "XY");
    TEST_ASSERT_NOT_NULL(aux);
    const char *aux_str = bam_aux2Z(aux);
    TEST_ASSERT_NOT_NULL(aux_str);
    TEST_ASSERT_EQUAL_STRING("foo", aux_str);

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fastq_record_to_bam_converts_long_quality_string)
{
    char name[] = "read1";
    char sequence[65];
    char quality[65];
    memset(sequence, 'A', sizeof(sequence) - 1);
    memset(quality, 'I', sizeof(quality) - 1);
    sequence[sizeof(sequence) - 1] = '\0';
    quality[sizeof(quality) - 1] = '\0';

    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_record_to_bam(bam,
                                 name, strlen(name), NULL, 0,
                                 sequence, quality, strlen(sequence)));
    uint8_t *qual = bam_get_qual(bam);
    for (size_t i = 0; i < strlen(quality); ++i) {
        TEST_ASSERT_EQUAL_UINT8(40, qual[i]);
    }

    bam_destroy1(bam);
}

FC_TEST(bam_utils_suite, test_fastq_record_to_bam_matches_htslib_for_iupac_sequence)
{
    char name[] = "iupac-read";
    char sequence[] = "ACGTNacgtnMRWSYKVHDBX=mrwsykvhdbxA";
    char quality[sizeof(sequence)];
    memset(quality, 'I', sizeof(quality) - 1);
    quality[sizeof(quality) - 1] = '\0';

    bam1_t *actual = bam_init1();
    bam1_t *expected = bam_init1();
    TEST_ASSERT_NOT_NULL(actual);
    TEST_ASSERT_NOT_NULL(expected);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_record_to_bam(actual,
                                 name, strlen(name), NULL, 0,
                                 sequence, quality, strlen(sequence)));
    TEST_ASSERT_TRUE(bam_set1(
        expected,
        strlen(name) + 1, name,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        strlen(sequence), sequence, quality,
        0) >= 0);
    uint8_t *expected_qual = bam_get_qual(expected);
    for (size_t i = 0; i < strlen(quality); ++i) {
        expected_qual[i] -= 33;
    }

    TEST_ASSERT_EQUAL_MEMORY(&expected->core, &actual->core, sizeof(actual->core));
    TEST_ASSERT_EQUAL_INT(expected->l_data, actual->l_data);
    TEST_ASSERT_EQUAL_MEMORY(expected->data, actual->data, (size_t)actual->l_data);

    bam_destroy1(actual);
    bam_destroy1(expected);
}

FC_TEST(bam_utils_suite, test_fastq_record_to_bam_matches_htslib_for_long_common_sequence)
{
    char name[] = "common-read";
    char sequence[] =
        "ACGTNacgtnACGTNacgtnACGTNacgtnAC"
        "GTNacgtnACGTNacgtnACGTNacgtnACGTN";
    char quality[sizeof(sequence)];
    memset(quality, 'I', sizeof(quality) - 1);
    quality[sizeof(quality) - 1] = '\0';

    bam1_t *actual = bam_init1();
    bam1_t *expected = bam_init1();
    TEST_ASSERT_NOT_NULL(actual);
    TEST_ASSERT_NOT_NULL(expected);
    TEST_ASSERT_EQUAL_INT(0, fc_fastq_record_to_bam(actual,
                                 name, strlen(name), NULL, 0,
                                 sequence, quality, strlen(sequence)));
    TEST_ASSERT_TRUE(bam_set1(
        expected,
        strlen(name) + 1, name,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        strlen(sequence), sequence, quality,
        0) >= 0);
    uint8_t *expected_qual = bam_get_qual(expected);
    for (size_t i = 0; i < strlen(quality); ++i) {
        expected_qual[i] -= 33;
    }

    TEST_ASSERT_EQUAL_MEMORY(&expected->core, &actual->core, sizeof(actual->core));
    TEST_ASSERT_EQUAL_INT(expected->l_data, actual->l_data);
    TEST_ASSERT_EQUAL_MEMORY(expected->data, actual->data, (size_t)actual->l_data);

    bam_destroy1(actual);
    bam_destroy1(expected);
}

FC_TEST(bam_utils_suite, test_fastq_record_to_bam_accepts_borrowed_fields)
{
    static const char input[] = "read42 comment=value\nACGTN\nIIIII";
    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);

    TEST_ASSERT_EQUAL_INT(0, fc_fastq_record_to_bam(bam, input, 6, input + 7, 13, input + 21, input + 27, 5));
    TEST_ASSERT_EQUAL_STRING("read42", bam_get_qname(bam));
    TEST_ASSERT_EQUAL_INT(5, bam->core.l_qseq);
    TEST_ASSERT_EQUAL_UINT8(40, bam_get_qual(bam)[0]);

    bam_destroy1(bam);
}

FC_TEST(bam_utils_suite, test_bam_input_info_retains_matching_reference_header)
{
    sam_hdr_t *first = header_from_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@SQ\tSN:chr2\tLN:200\n");
    sam_hdr_t *second = header_from_text(
        "@HD\tVN:1.6\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@SQ\tSN:chr2\tLN:200\n");
    fc_bam_input_info info = { 0 };

    TEST_ASSERT_EQUAL_INT(FC_BAM_INPUT_OK, fc_bam_input_info_add_header(&info, first));
    TEST_ASSERT_TRUE(info.has_references);
    TEST_ASSERT_NOT_NULL(info.reference_header);
    TEST_ASSERT_TRUE(fc_bam_header_is_coordinate_sorted(info.reference_header));
    TEST_ASSERT_EQUAL_INT(FC_BAM_INPUT_OK, fc_bam_input_info_add_header(&info, second));

    sam_hdr_destroy(first);
    sam_hdr_destroy(second);
    fc_bam_input_info_destroy(&info);
}

FC_TEST(bam_utils_suite, test_bam_input_info_rejects_different_reference_dictionary)
{
    sam_hdr_t *first = header_from_text("@SQ\tSN:chr1\tLN:100\n");
    sam_hdr_t *different_length = header_from_text("@SQ\tSN:chr1\tLN:101\n");
    sam_hdr_t *different_order = header_from_text(
        "@SQ\tSN:chr2\tLN:200\n"
        "@SQ\tSN:chr1\tLN:100\n");
    fc_bam_input_info info = { 0 };

    TEST_ASSERT_EQUAL_INT(FC_BAM_INPUT_OK, fc_bam_input_info_add_header(&info, first));
    TEST_ASSERT_EQUAL_INT(
        FC_BAM_INPUT_REFERENCE_MISMATCH,
        fc_bam_input_info_add_header(&info, different_length));
    TEST_ASSERT_EQUAL_INT(
        FC_BAM_INPUT_REFERENCE_MISMATCH,
        fc_bam_input_info_add_header(&info, different_order));
    TEST_ASSERT_FALSE(fc_bam_header_is_coordinate_sorted(info.reference_header));

    sam_hdr_destroy(first);
    sam_hdr_destroy(different_length);
    sam_hdr_destroy(different_order);
    fc_bam_input_info_destroy(&info);
}

FC_TEST(bam_utils_suite, test_bam_coordinate_order_rejects_backtracking_and_mapped_after_unmapped)
{
    bam1_t record = { 0 };
    fc_bam_coordinate_order order = { 0 };

    record.core.tid = 0;
    record.core.pos = 100;
    TEST_ASSERT_TRUE(fc_bam_coordinate_order_check(&order, &record));
    record.core.pos = 101;
    TEST_ASSERT_TRUE(fc_bam_coordinate_order_check(&order, &record));
    record.core.tid = 1;
    record.core.pos = 0;
    TEST_ASSERT_TRUE(fc_bam_coordinate_order_check(&order, &record));

    record.core.flag = BAM_FUNMAP;
    TEST_ASSERT_TRUE(fc_bam_coordinate_order_check(&order, &record));
    record.core.flag = 0;
    TEST_ASSERT_FALSE(fc_bam_coordinate_order_check(&order, &record));

    order = (fc_bam_coordinate_order){ 0 };
    record.core.tid = 1;
    record.core.pos = 100;
    TEST_ASSERT_TRUE(fc_bam_coordinate_order_check(&order, &record));
    record.core.tid = 0;
    TEST_ASSERT_FALSE(fc_bam_coordinate_order_check(&order, &record));
}

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_sam_comment_consistency)
{
    const char *runid = "12345678-1234-1234-1234-123456789012";
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr,
        "RG:Z:12345678-1234-1234-1234-123456789012_sup@v1_barcode07\t"
        "SM:Z:barcode07\tpi:Z:parent-uuid\tML:B:C,1,2\t"
        "st:Z:2023-01-01T00:00:00Z\tdu:f:1.25\trn:i:7\tch:i:42\t"
        "PU:Z:FLO-MIN114\tLB:Z:sample123\tal:Z:alias1\tDS:Z:gpu:gpu0");

    uint8_t *aux_rg = bam_aux_get(fr->bam, "RG");
    TEST_ASSERT_NOT_NULL(aux_rg);
    TEST_ASSERT_EQUAL_STRING("12345678-1234-1234-1234-123456789012_sup@v1_barcode07", bam_aux2Z(aux_rg));

    uint8_t *aux_sm = bam_aux_get(fr->bam, "SM");
    TEST_ASSERT_NOT_NULL(aux_sm);
    TEST_ASSERT_EQUAL_STRING("barcode07", bam_aux2Z(aux_sm));

    uint8_t *aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("parent-uuid", bam_aux2Z(aux_pi));

    uint8_t *aux_ml = bam_aux_get(fr->bam, "ML");
    TEST_ASSERT_NOT_NULL(aux_ml);
    TEST_ASSERT_EQUAL_UINT32(2, bam_auxB_len(aux_ml));
    TEST_ASSERT_EQUAL_INT64(1, bam_auxB2i(aux_ml, 0));
    TEST_ASSERT_EQUAL_INT64(2, bam_auxB2i(aux_ml, 1));

    uint8_t *aux_st = bam_aux_get(fr->bam, "st");
    TEST_ASSERT_NOT_NULL(aux_st);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", bam_aux2Z(aux_st));

    uint8_t *aux_du = bam_aux_get(fr->bam, "du");
    TEST_ASSERT_NOT_NULL(aux_du);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, bam_aux2f(aux_du));

    uint8_t *aux_rn = bam_aux_get(fr->bam, "rn");
    TEST_ASSERT_NOT_NULL(aux_rn);
    TEST_ASSERT_EQUAL_INT(7, bam_aux2i(aux_rn));

    uint8_t *aux_ch = bam_aux_get(fr->bam, "ch");
    TEST_ASSERT_NOT_NULL(aux_ch);
    TEST_ASSERT_EQUAL_INT(42, bam_aux2i(aux_ch));

    uint8_t *aux_pu = bam_aux_get(fr->bam, "PU");
    TEST_ASSERT_NOT_NULL(aux_pu);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", bam_aux2Z(aux_pu));

    uint8_t *aux_lb = bam_aux_get(fr->bam, "LB");
    TEST_ASSERT_NOT_NULL(aux_lb);
    TEST_ASSERT_EQUAL_STRING("sample123", bam_aux2Z(aux_lb));

    uint8_t *aux_al = bam_aux_get(fr->bam, "al");
    TEST_ASSERT_NOT_NULL(aux_al);
    TEST_ASSERT_EQUAL_STRING("alias1", bam_aux2Z(aux_al));

    uint8_t *aux_ds = bam_aux_get(fr->bam, "DS");
    TEST_ASSERT_NOT_NULL(aux_ds);
    TEST_ASSERT_EQUAL_STRING("gpu:gpu0", bam_aux2Z(aux_ds));

    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING("sup@v1", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("", fr->meta->modcallers);
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

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_key_value_fallback_consistency)
{
    const char *runid = "e4994c62-93f9-439a-bc8f-d20c95a137a5";
    const char *model = "rna004_130bps_fast@v5.1.0";
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr,
        "runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 "
        "basecall_model_version_id=rna004_130bps_fast@v5.1.0 "
        "flow_cell_id=FLO-MIN114 sample_id=sample123 "
        "barcode=barcode01 barcode_alias=alias1 read=7 ch=42 "
        "start_time=2023-01-01T00:00:00Z du=2.5 parent_read_id=pr-001 "
        "basecall_gpu=gpu0 empty= unknown=foo");

    uint8_t *aux_rg = bam_aux_get(fr->bam, "RG");
    TEST_ASSERT_NOT_NULL(aux_rg);
    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode01",
        bam_aux2Z(aux_rg));

    uint8_t *aux_pu = bam_aux_get(fr->bam, "PU");
    TEST_ASSERT_NOT_NULL(aux_pu);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", bam_aux2Z(aux_pu));

    uint8_t *aux_lb = bam_aux_get(fr->bam, "LB");
    TEST_ASSERT_NOT_NULL(aux_lb);
    TEST_ASSERT_EQUAL_STRING("sample123", bam_aux2Z(aux_lb));

    uint8_t *aux_al = bam_aux_get(fr->bam, "al");
    TEST_ASSERT_NOT_NULL(aux_al);
    TEST_ASSERT_EQUAL_STRING("alias1", bam_aux2Z(aux_al));

    uint8_t *aux_pi = bam_aux_get(fr->bam, "pi");
    TEST_ASSERT_NOT_NULL(aux_pi);
    TEST_ASSERT_EQUAL_STRING("pr-001", bam_aux2Z(aux_pi));

    uint8_t *aux_sm = bam_aux_get(fr->bam, "SM");
    TEST_ASSERT_NOT_NULL(aux_sm);
    TEST_ASSERT_EQUAL_STRING("barcode01", bam_aux2Z(aux_sm));

    uint8_t *aux_rn = bam_aux_get(fr->bam, "rn");
    TEST_ASSERT_NOT_NULL(aux_rn);
    TEST_ASSERT_EQUAL_INT(7, bam_aux2i(aux_rn));

    uint8_t *aux_ch = bam_aux_get(fr->bam, "ch");
    TEST_ASSERT_NOT_NULL(aux_ch);
    TEST_ASSERT_EQUAL_INT(42, bam_aux2i(aux_ch));

    uint8_t *aux_st = bam_aux_get(fr->bam, "st");
    TEST_ASSERT_NOT_NULL(aux_st);
    TEST_ASSERT_EQUAL_STRING("2023-01-01T00:00:00Z", bam_aux2Z(aux_st));

    uint8_t *aux_du = bam_aux_get(fr->bam, "du");
    TEST_ASSERT_NOT_NULL(aux_du);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 2.5, bam_aux2f(aux_du));

    uint8_t *aux_ds = bam_aux_get(fr->bam, "DS");
    TEST_ASSERT_NOT_NULL(aux_ds);
    TEST_ASSERT_EQUAL_STRING("gpu:gpu0", bam_aux2Z(aux_ds));

    uint8_t *aux_co = bam_aux_get(fr->bam, "CO");
    TEST_ASSERT_NOT_NULL(aux_co);
    TEST_ASSERT_EQUAL_STRING("empty= unknown=foo", bam_aux2Z(aux_co));

    TEST_ASSERT_EQUAL_STRING(runid, fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING(model, fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("", fr->meta->modcallers);
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

FC_TEST(bam_utils_suite, test_fc_read_prefers_rg_ds_models_from_header)
{
    bam_hdr_t *hdr = NULL;
    fc_bam_rg_cache *cache = NULL;
    fc_read *fr = read_from_rg_header_ds(
        "rg1",
        "basecall_model=rna004_130bps_sup@2025-06-12_n1mpU-integrity "
        "modbase_models=rna004_130bps_hac@v5.1.0_inosine_m6A@v1,"
        "rna004_130bps_hac@v5.2.0_pseU@v1",
        &hdr,
        &cache);

    TEST_ASSERT_EQUAL_STRING("rna004_130bps_sup@2025-06-12_n1mpU-integrity", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING(
        "rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1",
        fr->meta->modcallers);

    destroy_read(fr);
    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_utils_suite, test_fc_read_parses_rg_ds_basecall_only)
{
    bam_hdr_t *hdr = NULL;
    fc_bam_rg_cache *cache = NULL;
    fc_read *fr = read_from_rg_header_ds(
        "rg1",
        "basecall_model=rna004_130bps_fast@v5.1.0",
        &hdr,
        &cache);

    TEST_ASSERT_EQUAL_STRING("rna004_130bps_fast@v5.1.0", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("", fr->meta->modcallers);
    destroy_read(fr);
    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_utils_suite, test_fc_read_parses_rg_ds_single_modbase_model)
{
    bam_hdr_t *hdr = NULL;
    fc_bam_rg_cache *cache = NULL;
    fc_read *fr = read_from_rg_header_ds(
        "rg1",
        "basecall_model=rna004_130bps_sup@2024-06-05_canonical "
        "modbase_models=modbase_model_version_id",
        &hdr,
        &cache);

    TEST_ASSERT_EQUAL_STRING("rna004_130bps_sup@2024-06-05_canonical", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("modbase_model_version_id", fr->meta->modcallers);
    destroy_read(fr);
    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_utils_suite, test_fc_read_parses_rg_ds_multiple_modbase_models_as_one_string)
{
    bam_hdr_t *hdr = NULL;
    fc_bam_rg_cache *cache = NULL;
    fc_read *fr = read_from_rg_header_ds(
        "rg1",
        "basecall_model=rna004_130bps_fast@v5.1.0 "
        "modbase_models=rna004_130bps_hac@v5.1.0_inosine_m6A@v1,"
        "rna004_130bps_hac@v5.2.0_pseU@v1",
        &hdr,
        &cache);

    TEST_ASSERT_EQUAL_STRING("rna004_130bps_fast@v5.1.0", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING(
        "rna004_130bps_hac@v5.1.0_inosine_m6A@v1,rna004_130bps_hac@v5.2.0_pseU@v1",
        fr->meta->modcallers);
    destroy_read(fr);
    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_utils_suite, test_fc_read_parses_rg_ds_integrity_suffix_variants)
{
    bam_hdr_t *canonical_hdr = NULL;
    fc_bam_rg_cache *canonical_cache = NULL;
    fc_read *canonical = read_from_rg_header_ds(
        "rg1",
        "basecall_model=rna004_130bps_sup@2025-06-12_canonical-integrity "
        "modbase_models=modbase_model_version_id",
        &canonical_hdr,
        &canonical_cache);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_sup@2025-06-12_canonical-integrity", canonical->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("modbase_model_version_id", canonical->meta->modcallers);
    destroy_read(canonical);
    fc_bam_rg_cache_destroy(canonical_cache);
    sam_hdr_destroy(canonical_hdr);

    bam_hdr_t *n1mpu_hdr = NULL;
    fc_bam_rg_cache *n1mpu_cache = NULL;
    fc_read *n1mpu = read_from_rg_header_ds(
        "rg2",
        "basecall_model=rna004_130bps_sup@2025-06-12_n1mpU-integrity "
        "modbase_models=modbase_model_version_id",
        &n1mpu_hdr,
        &n1mpu_cache);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_sup@2025-06-12_n1mpU-integrity", n1mpu->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("modbase_model_version_id", n1mpu->meta->modcallers);
    destroy_read(n1mpu);
    fc_bam_rg_cache_destroy(n1mpu_cache);
    sam_hdr_destroy(n1mpu_hdr);
}

FC_TEST(bam_utils_suite, test_bam_rg_cache_lookup_and_rg_field_fallbacks)
{
    const char *hdr_text =
        "@HD\tVN:1.6\n"
        "@RG\tID:rg1\tDT:2025-01-06T10:06:36Z\tLB:sampleA\tPU:FLO-MIN114\tPM:MIN-001\tSM:barcode07\tal:alias07\tDS:runid=e4994c62-93f9-439a-bc8f-d20c95a137a5 basecall_model=rna004_130bps_fast@v5.1.0 modbase_models=rna004_130bps_hac@v5.1.0_inosine_m6A@v1\n"
        "@RG\tID:rg2\tDS:basecall_model=rna004_130bps_sup@2025-06-12_canonical-integrity\n";
    bam_hdr_t *hdr = sam_hdr_parse((int)strlen(hdr_text), hdr_text);
    TEST_ASSERT_NOT_NULL(hdr);

    fc_bam_rg_cache *cache = fc_bam_rg_cache_create(hdr);
    TEST_ASSERT_NOT_NULL(cache);

    const fc_bam_rg_entry *rg1 = fc_bam_rg_cache_find(cache, "rg1");
    TEST_ASSERT_NOT_NULL(rg1);
    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5", rg1->runid);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_fast@v5.1.0", rg1->basecall_model);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_hac@v5.1.0_inosine_m6A@v1", rg1->modbase_models);
    TEST_ASSERT_EQUAL_STRING("2025-01-06T10:06:36Z", rg1->dt);
    TEST_ASSERT_EQUAL_STRING("sampleA", rg1->lb);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", rg1->pu);
    TEST_ASSERT_EQUAL_STRING("MIN-001", rg1->pm);
    TEST_ASSERT_EQUAL_STRING("barcode07", rg1->sm);
    TEST_ASSERT_EQUAL_STRING("alias07", rg1->al);
    TEST_ASSERT_NULL(fc_bam_rg_cache_find(cache, "missing"));

    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    const char *name = "read1";
    const char *seq = "ACGT";
    const char *qual_ascii = "####";
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
    TEST_ASSERT_TRUE(bam_aux_update_str(fr->bam, "RG", 4, "rg1") >= 0);
    fc_bam_ensure_mean_quality(fr->bam, false);
    fc_read_reset_from_bam(fr, hdr, cache, NULL);

    TEST_ASSERT_EQUAL_STRING("e4994c62-93f9-439a-bc8f-d20c95a137a5", fr->meta->runid);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_fast@v5.1.0", fr->meta->basecaller);
    TEST_ASSERT_EQUAL_STRING("rna004_130bps_hac@v5.1.0_inosine_m6A@v1", fr->meta->modcallers);
    TEST_ASSERT_EQUAL_STRING("2025-01-06T10:06:36Z", fr->meta->start_time);
    TEST_ASSERT_EQUAL_STRING("sampleA", fr->meta->sample_id);
    TEST_ASSERT_EQUAL_STRING("FLO-MIN114", fr->meta->flow_cell_id);
    TEST_ASSERT_EQUAL_STRING("MIN-001", fr->meta->device_id);
    TEST_ASSERT_EQUAL_STRING("barcode07", fr->meta->barcode);
    TEST_ASSERT_EQUAL_size_t(7, fr->meta->ibarcode);
    TEST_ASSERT_EQUAL_STRING("alias07", fr->meta->barcode_alias);

    destroy_read(fr);
    fc_bam_rg_cache_destroy(cache);
    sam_hdr_destroy(hdr);
}

FC_TEST(bam_utils_suite, test_fc_read_preserves_original_comment)
{
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "raw	comment  with  spaces");
    TEST_ASSERT_NOT_NULL(fr->meta->orig_comment);
    TEST_ASSERT_EQUAL_STRING("raw	comment  with  spaces", fr->meta->orig_comment);

    destroy_read(fr);
}

FC_TEST(bam_utils_suite, test_fc_read_reset_from_kseq_reuses_preallocated_bam)
{
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    read_from_comment(fr, "XY:Z:foo");
    bam1_t *first_bam = fr->bam;
    TEST_ASSERT_NOT_NULL(first_bam);

    read_from_comment(fr, "XY:Z:bar");
    TEST_ASSERT_EQUAL_PTR(first_bam, fr->bam);

    uint8_t *aux = bam_aux_get(fr->bam, "XY");
    TEST_ASSERT_NOT_NULL(aux);
    TEST_ASSERT_EQUAL_STRING("bar", bam_aux2Z(aux));
    destroy_read(fr);
}
