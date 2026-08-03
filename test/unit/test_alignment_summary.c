#include <string.h>

#include "core/fc_alignment_summary.h"
#include "core/fc_writer.h"
#include "main.h"
#include "unity.h"

static bam1_t alignment_test_bam(int32_t tid, uint16_t flag)
{
    bam1_t bam = { 0 };
    bam.core.tid = tid;
    bam.core.flag = flag;
    return bam;
}

FC_TEST(alignment_summary_suite, test_alignment_summary_counts_all_flag_categories)
{
    fc_alignment_summary *summary = fc_alignment_summary_create(2, false, false, 1);
    TEST_ASSERT_NOT_NULL(summary);

    bam1_t primary = alignment_test_bam(0, 0);
    bam1_t secondary = alignment_test_bam(0, BAM_FSECONDARY | BAM_FDUP);
    bam1_t supplementary = alignment_test_bam(1, BAM_FSUPPLEMENTARY | BAM_FQCFAIL);
    bam1_t unmapped = alignment_test_bam(-1, BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP);
    bam1_t secondary_unmapped = alignment_test_bam(-1, BAM_FUNMAP | BAM_FSECONDARY);
    fc_alignment_summary_add_read(summary, 0, &primary, 1);
    fc_alignment_summary_add_read(summary, 0, &secondary, 0);
    fc_alignment_summary_add_read(summary, 0, &supplementary, -1);
    fc_alignment_summary_add_read(summary, 0, &unmapped, 0);
    fc_alignment_summary_add_read(summary, 0, &secondary_unmapped, 0);

    const fc_alignment_counts *counts = fc_alignment_summary_global(summary);
    TEST_ASSERT_NOT_NULL(counts);
    TEST_ASSERT_TRUE(counts->seen);
    TEST_ASSERT_EQUAL_UINT64(2, counts->reference_counts[0][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(1, counts->reference_counts[0][FC_FLAGSTAT_PRIMARY]);
    TEST_ASSERT_EQUAL_UINT64(1, counts->reference_counts[0][FC_FLAGSTAT_SECONDARY]);
    TEST_ASSERT_EQUAL_UINT64(1, counts->reference_counts[1][FC_FLAGSTAT_SUPPLEMENTARY]);
    TEST_ASSERT_EQUAL_UINT64(2, counts->unmapped_counts[FC_FLAGSTAT_UNMAPPED]);

    fc_alignment_totals totals = { 0 };
    fc_alignment_counts_materialize(counts, &totals);
    TEST_ASSERT_EQUAL_UINT64(5, totals.record_counts[FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(1, totals.record_counts[FC_FLAGSTAT_PRIMARY]);
    TEST_ASSERT_EQUAL_UINT64(2, totals.record_counts[FC_FLAGSTAT_SECONDARY]);
    TEST_ASSERT_EQUAL_UINT64(1, totals.record_counts[FC_FLAGSTAT_SUPPLEMENTARY]);
    TEST_ASSERT_EQUAL_UINT64(2, totals.record_counts[FC_FLAGSTAT_UNMAPPED]);
    TEST_ASSERT_EQUAL_UINT64(2, totals.record_counts[FC_FLAGSTAT_QCFAIL]);
    TEST_ASSERT_EQUAL_UINT64(2, totals.record_counts[FC_FLAGSTAT_DUPLICATE]);
    TEST_ASSERT_EQUAL_UINT64(1, totals.record_counts[FC_FLAGSTAT_DUPLEX]);
    TEST_ASSERT_EQUAL_UINT64(1, totals.record_counts[FC_FLAGSTAT_DUPLEX_FORMING]);
    TEST_ASSERT_EQUAL_UINT64(3, totals.mapped_record_count);
    TEST_ASSERT_EQUAL_UINT64(1, totals.primary_mapped_record_count);
    TEST_ASSERT_EQUAL_UINT64(1, totals.mapped_read_count);
    TEST_ASSERT_EQUAL_UINT64(2, totals.read_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, totals.mapping_rate);

    fc_alignment_summary_destroy(summary);
}

FC_TEST(alignment_summary_suite, test_alignment_summary_tracks_global_and_barcode_counts)
{
    fc_alignment_summary *summary = fc_alignment_summary_create(2, false, true, 3);
    TEST_ASSERT_NOT_NULL(summary);

    bam1_t ref0 = alignment_test_bam(0, 0);
    bam1_t ref1 = alignment_test_bam(1, 0);
    bam1_t unmapped = alignment_test_bam(-1, BAM_FUNMAP);
    bam1_t secondary = alignment_test_bam(1, BAM_FSECONDARY);
    fc_alignment_summary_add_read(summary, 0, &ref0, 0);
    fc_alignment_summary_add_read(summary, 2, &ref1, 0);
    fc_alignment_summary_add_read(summary, 2, &unmapped, 0);
    fc_alignment_summary_add_read(summary, 9, &secondary, 0);

    const fc_alignment_counts *barcode0 = fc_alignment_summary_barcode(summary, 0);
    const fc_alignment_counts *barcode2 = fc_alignment_summary_barcode(summary, 2);
    TEST_ASSERT_NOT_NULL(barcode0);
    TEST_ASSERT_NOT_NULL(barcode2);
    TEST_ASSERT_NULL(fc_alignment_summary_barcode(summary, 1));
    TEST_ASSERT_EQUAL_PTR(barcode0, fc_alignment_summary_for_output(summary, 0));
    TEST_ASSERT_EQUAL_UINT64(1, barcode0->reference_counts[0][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(1, barcode2->reference_counts[1][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(1, barcode2->unmapped_counts[FC_FLAGSTAT_TOTAL]);

    fc_alignment_totals totals = { 0 };
    fc_alignment_counts_materialize(fc_alignment_summary_global(summary), &totals);
    TEST_ASSERT_EQUAL_UINT64(4, totals.record_counts[FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(3, totals.mapped_record_count);
    TEST_ASSERT_EQUAL_UINT64(2, totals.primary_mapped_record_count);
    TEST_ASSERT_EQUAL_UINT64(3, totals.read_count);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0 / 3.0, totals.mapping_rate);

    fc_alignment_summary_destroy(summary);
}

FC_TEST(alignment_summary_suite, test_alignment_summary_region_mode_assigns_global_and_barcode)
{
    fc_alignment_summary *summary = fc_alignment_summary_create(3, true, true, 2);
    TEST_ASSERT_NOT_NULL(summary);

    bam1_t mapped = alignment_test_bam(77, 0);
    fc_alignment_summary_add_read(summary, 1, &mapped, 0);

    const fc_alignment_counts *global = fc_alignment_summary_global(summary);
    const fc_alignment_counts *barcode = fc_alignment_summary_barcode(summary, 1);
    TEST_ASSERT_EQUAL_UINT64(1, global->n_refs);
    TEST_ASSERT_EQUAL_UINT64(1, global->reference_counts[0][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_NOT_NULL(barcode);
    TEST_ASSERT_EQUAL_UINT64(1, barcode->n_refs);
    TEST_ASSERT_EQUAL_UINT64(1, barcode->reference_counts[0][FC_FLAGSTAT_TOTAL]);

    fc_alignment_summary_destroy(summary);
}

FC_TEST(alignment_summary_suite, test_alignment_summary_ignores_unknown_reference)
{
    fc_alignment_summary *summary = fc_alignment_summary_create(1, false, false, 1);
    TEST_ASSERT_NOT_NULL(summary);

    bam1_t mapped = alignment_test_bam(7, 0);
    fc_alignment_summary_add_read(summary, 0, &mapped, 0);

    const fc_alignment_counts *counts = fc_alignment_summary_global(summary);
    TEST_ASSERT_FALSE(counts->seen);
    fc_alignment_totals totals = { 0 };
    fc_alignment_counts_materialize(counts, &totals);
    TEST_ASSERT_EQUAL_UINT64(0, totals.record_counts[FC_FLAGSTAT_TOTAL]);

    fc_alignment_summary_destroy(summary);
}

FC_TEST(alignment_summary_suite, test_writer_updates_shared_alignment_summary_once)
{
    fc_alignment_summary *summary = fc_alignment_summary_create(1, false, true, 2);
    TEST_ASSERT_NOT_NULL(summary);
    fc_writer_config config = {
        .use_barcodes = true,
        .alignment_summary = summary
    };
    fc_writer *writer = fc_writer_create(&config);
    TEST_ASSERT_NOT_NULL(writer);

    bam1_t bam = alignment_test_bam(0, 0);
    read_meta meta = {
        .ibarcode = 1,
        .duplex_code = 1
    };
    fc_read read = {
        .bam = &bam,
        .meta = &meta
    };
    fc_writer_add_read(writer, &read, "input.bam", NULL);

    const fc_alignment_counts *global = fc_alignment_summary_global(summary);
    const fc_alignment_counts *barcode = fc_alignment_summary_barcode(summary, 1);
    TEST_ASSERT_EQUAL_UINT64(1, global->reference_counts[0][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_NOT_NULL(barcode);
    TEST_ASSERT_EQUAL_UINT64(1, barcode->reference_counts[0][FC_FLAGSTAT_TOTAL]);
    TEST_ASSERT_EQUAL_UINT64(1, barcode->reference_counts[0][FC_FLAGSTAT_DUPLEX]);

    fc_writer_destroy(writer);
    fc_alignment_summary_destroy(summary);
}
