#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "main.h"
#include "htslib/sam.h"

#include "core/fc_stats_engine.h"
#include "core/fc_bam_utils.h"

static size_t sum_counts(const read_stats* stats) {
    size_t total = 0;
    if (stats == NULL) return 0;
    for (size_t i = 0; i < stats->n; ++i) {
        total += stats->counts[i];
    }
    return total;
}

static bam1_t* build_unmapped_read(const char* qname, const char* seq, const char* qual) {
    int32_t l_seq = (int32_t)strlen(seq);
    int32_t l_qname = (int32_t)strlen(qname) + 1;

    bam1_t* b = bam_init1();
    bam_set1(
        b, l_qname, qname,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        l_seq, seq, qual,
        0);

    uint8_t* q = bam_get_qual(b);
    for (int i = 0; i < l_seq; ++i) {
        q[i] = (uint8_t)(qual[i] - 33);
    }

    return b;
}

static bam1_t* build_mapped_read(const char* qname, const char* seq, const char* qual, int32_t tid, int32_t pos) {
    int32_t l_seq = (int32_t)strlen(seq);
    int32_t l_qname = (int32_t)strlen(qname) + 1;
    uint32_t cigar = bam_cigar_gen(l_seq, BAM_CMATCH);

    bam1_t* b = bam_init1();
    bam_set1(
        b, l_qname, qname,
        0, tid, pos, 60,
        1, &cigar,
        -1, -1, 0,
        l_seq, seq, qual,
        0);

    uint8_t* q = bam_get_qual(b);
    for (int i = 0; i < l_seq; ++i) {
        q[i] = (uint8_t)(qual[i] - 33);
    }
    TEST_ASSERT_TRUE(bam_aux_update_int(b, "NM", 0) >= 0);

    return b;
}

FC_TEST(stats_engine_suite, test_stats_engine_unmapped_tracked_in_canonical_stats)
{
    bam1_t* b = build_unmapped_read("r1", "ACGT", "!!!!");
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    bam_destroy1(fr->bam);
    fr->bam = b;
    fc_bam_ensure_mean_quality(b, false);
    fc_read_reset_from_bam(fr, NULL, fc_test_rg_cache(), NULL);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 0.0f, fr->meta->mean_q);

    fc_stats_engine* engine = fc_stats_engine_create(false, false, 0.0f, 0.0f, false);
    TEST_ASSERT_EQUAL_INT(FC_STATS_UNMAPPED, fc_stats_engine_add_read(engine, fr, NULL, NULL));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->qual_stats));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(fr);
}

FC_TEST(stats_engine_suite, test_stats_engine_aligned_metrics)
{
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "100", NULL);

    bam1_t* b = build_mapped_read("r1", "ACGTACGTAC", "!!!!!!!!!!", 0, 5);
    fc_read *fr = fc_read_init();
    TEST_ASSERT_NOT_NULL(fr);
    bam_destroy1(fr->bam);
    fr->bam = b;
    fc_bam_ensure_mean_quality(b, false);
    fc_read_reset_from_bam(fr, hdr, fc_test_rg_cache(), NULL);
    TEST_ASSERT_FLOAT_WITHIN(1e-3, 0.0f, fr->meta->mean_q);

    fc_stats_engine* engine = fc_stats_engine_create(true, false, 0.0f, 0.0f, false);
    fc_stats_aligned metrics = {0};
    TEST_ASSERT_EQUAL_INT(FC_STATS_ALIGNED, fc_stats_engine_add_read(engine, fr, hdr, &metrics));

    TEST_ASSERT_EQUAL_UINT32(10, metrics.read_length);
    TEST_ASSERT_EQUAL_UINT32(5, metrics.rstart);
    TEST_ASSERT_EQUAL_UINT32(15, metrics.rend);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.aligned_ref_len);
    TEST_ASSERT_EQUAL_UINT32(100, metrics.ref_length);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.qstart);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.qend);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.match);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.ins);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.delt);
    TEST_ASSERT_EQUAL_UINT32(0, metrics.sub);
    TEST_ASSERT_EQUAL_UINT32(10, metrics.length);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.coverage);
    TEST_ASSERT_EQUAL_FLOAT(10.0f, metrics.ref_cover);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.iden);
    TEST_ASSERT_EQUAL_FLOAT(100.0f, metrics.acc);
    TEST_ASSERT_EQUAL_CHAR('+', metrics.direction);

    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->qual_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->acc_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->cov_stats));
    TEST_ASSERT_EQUAL_UINT64(1, engine->acc_stats->counts[engine->acc_stats->n - 1]);
    TEST_ASSERT_EQUAL_UINT64(1, engine->cov_stats->counts[engine->cov_stats->n - 1]);

    fc_stats_engine_destroy(engine);
    fc_read_destroy(fr);
    sam_hdr_destroy(hdr);
}

FC_TEST(stats_engine_suite, test_stats_engine_alignment_histogram_populations)
{
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "100", NULL);

    bam1_t* mapped = build_mapped_read("mapped", "ACGTACGTAC", "!!!!!!!!!!", 0, 5);
    bam1_t* unmapped = build_unmapped_read("unmapped", "ACGT", "!!!!");
    fc_read* mapped_read = fc_read_init();
    fc_read* unmapped_read = fc_read_init();
    TEST_ASSERT_NOT_NULL(mapped_read);
    TEST_ASSERT_NOT_NULL(unmapped_read);
    bam_destroy1(mapped_read->bam);
    bam_destroy1(unmapped_read->bam);
    mapped_read->bam = mapped;
    unmapped_read->bam = unmapped;
    fc_bam_ensure_mean_quality(mapped, false);
    fc_bam_ensure_mean_quality(unmapped, false);
    fc_read_reset_from_bam(mapped_read, hdr, fc_test_rg_cache(), NULL);
    fc_read_reset_from_bam(unmapped_read, NULL, fc_test_rg_cache(), NULL);
    mapped_read->meta->poly_tail_length = 7;

    fc_stats_engine* engine = fc_stats_engine_create(true, true, 0.0f, 0.0f, false);
    TEST_ASSERT_EQUAL_INT(FC_STATS_ALIGNED, fc_stats_engine_add_read(engine, mapped_read, hdr, NULL));
    TEST_ASSERT_EQUAL_INT(FC_STATS_UNMAPPED, fc_stats_engine_add_read(engine, unmapped_read, NULL, NULL));

    TEST_ASSERT_EQUAL_UINT64(1, engine->length_stats->counts[4]);
    TEST_ASSERT_EQUAL_UINT64(2, engine->qual_stats->counts[0]);
    TEST_ASSERT_EQUAL_UINT64(1, engine->polya_stats->counts[7]);
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->acc_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->cov_stats));
    TEST_ASSERT_EQUAL_UINT64(2, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(2, sum_counts(engine->qual_stats));
    TEST_ASSERT_EQUAL_UINT64(1, sum_counts(engine->polya_stats));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(mapped_read);
    fc_read_destroy(unmapped_read);
    sam_hdr_destroy(hdr);
}

FC_TEST(stats_engine_suite, test_stats_engine_skips_additional_alignments)
{
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "100", NULL);

    bam1_t* secondary = build_mapped_read("secondary", "ACGT", "IIII", 0, 5);
    bam1_t* supplementary = build_mapped_read("supplementary", "ACGT", "IIII", 0, 10);
    secondary->core.flag |= BAM_FSECONDARY;
    supplementary->core.flag |= BAM_FSUPPLEMENTARY;

    fc_read* secondary_read = fc_read_init();
    fc_read* supplementary_read = fc_read_init();
    TEST_ASSERT_NOT_NULL(secondary_read);
    TEST_ASSERT_NOT_NULL(supplementary_read);
    bam_destroy1(secondary_read->bam);
    bam_destroy1(supplementary_read->bam);
    secondary_read->bam = secondary;
    supplementary_read->bam = supplementary;
    fc_bam_ensure_mean_quality(secondary, false);
    fc_bam_ensure_mean_quality(supplementary, false);
    fc_read_reset_from_bam(secondary_read, hdr, fc_test_rg_cache(), NULL);
    fc_read_reset_from_bam(supplementary_read, hdr, fc_test_rg_cache(), NULL);

    fc_stats_engine* engine = fc_stats_engine_create(true, false, 0.0f, 0.0f, false);
    TEST_ASSERT_EQUAL_INT(FC_STATS_SKIP, fc_stats_engine_add_read(engine, secondary_read, hdr, NULL));
    TEST_ASSERT_EQUAL_INT(FC_STATS_SKIP, fc_stats_engine_add_read(engine, supplementary_read, hdr, NULL));
    TEST_ASSERT_EQUAL_UINT64(0, sum_counts(engine->length_stats));
    TEST_ASSERT_EQUAL_UINT64(0, sum_counts(engine->acc_stats));

    fc_stats_engine_destroy(engine);
    fc_read_destroy(secondary_read);
    fc_read_destroy(supplementary_read);
    sam_hdr_destroy(hdr);
}
