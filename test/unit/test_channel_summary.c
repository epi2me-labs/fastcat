#include <string.h>

#include "unity.h"
#include "main.h"
#include "core/fc_channel_summary.h"
#include "core/fc_read.h"

static fc_read *make_channel_test_read(const char *name)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);

    const char *seq = "ACGT";
    const char *qual_ascii = "IIII";
    TEST_ASSERT_TRUE(bam_set1(
        read->bam,
        strlen(name) + 1,
        name,
        BAM_FUNMAP,
        -1,
        -1,
        0,
        0,
        NULL,
        -1,
        -1,
        0,
        strlen(seq),
        seq,
        qual_ascii,
        0) >= 0);
    return read;
}

static void add_channel_test_read(
    fc_channel_summary *model,
    size_t barcode,
    size_t channel,
    size_t mux,
    size_t length,
    float mean_quality)
{
    bam1_t bam = { 0 };
    read_meta meta = {
        .channel = channel,
        .mux = mux,
        .mean_q = mean_quality
    };
    fc_read read = {
        .bam = &bam,
        .meta = &meta
    };
    bam.core.l_qseq = (int32_t)length;
    fc_channel_summary_add_read(model, barcode, &read);
}

FC_TEST(channel_summary_suite, test_channel_summary_accumulates_counts)
{
    fc_channel_summary *model = fc_channel_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_channel_test_read(model, 0, 10, 2, 100, 10.0f);
    add_channel_test_read(model, 0, 10, 2, 150, 20.0f);
    add_channel_test_read(model, 0, 10, 4, 75, 30.0f);

    fc_channel_summary_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_channel_summary_materialize_sparse(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(2, table.n);
    TEST_ASSERT_EQUAL_size_t(10, table.highest_channel);

    TEST_ASSERT_EQUAL_size_t(10, table.rows[0].channel);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].mux);
    TEST_ASSERT_EQUAL_size_t(2, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(250, table.rows[0].base_count);
    TEST_ASSERT_EQUAL_DOUBLE(15.0, table.rows[0].mean_quality);

    TEST_ASSERT_EQUAL_size_t(10, table.rows[1].channel);
    TEST_ASSERT_EQUAL_size_t(4, table.rows[1].mux);
    TEST_ASSERT_EQUAL_size_t(1, table.rows[1].read_count);
    TEST_ASSERT_EQUAL_size_t(75, table.rows[1].base_count);
    TEST_ASSERT_EQUAL_DOUBLE(30.0, table.rows[1].mean_quality);

    fc_channel_summary_table_destroy(&table);
    fc_channel_summary_destroy(model);
}

FC_TEST(channel_summary_suite, test_channel_summary_accumulates_alignment_concordance)
{
    fc_channel_summary *model = fc_channel_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    uint32_t cigar = bam_cigar_gen(10, BAM_CMATCH);
    bam1_t *bam = bam_init1();
    TEST_ASSERT_NOT_NULL(bam);
    TEST_ASSERT_TRUE(bam_set1(
        bam, 11, "aligned-qc", 0, 0, 0, 60, 1, &cigar, -1, -1, 0, 10, "ACGTACGTAC",
        "IIIIIIIIII", 0) >= 0);
    TEST_ASSERT_TRUE(bam_aux_update_int(bam, "NM", 1) >= 0);
    read_meta meta = {
        .channel = 10,
        .mux = 2,
        .mean_q = 20.0f,
        .edit_distance = 1
    };
    fc_read read = {
        .bam = bam,
        .meta = &meta
    };
    fc_channel_summary_add_read(model, 0, &read);

    fc_channel_summary_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_channel_summary_materialize_sparse(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(1, table.n);
    TEST_ASSERT_EQUAL_size_t(1, table.rows[0].concordance_read_count);
    TEST_ASSERT_EQUAL_DOUBLE(10.0, table.rows[0].mean_concordance_qscore);

    fc_channel_summary_table_destroy(&table);
    bam_destroy1(bam);
    fc_channel_summary_destroy(model);
}

FC_TEST(channel_summary_suite, test_channel_summary_barcode_state_is_isolated)
{
    fc_channel_summary *model = fc_channel_summary_create(3);
    TEST_ASSERT_NOT_NULL(model);

    add_channel_test_read(model, 1, 5, 1, 10, 10.0f);
    add_channel_test_read(model, 2, 7, 4, 20, 20.0f);

    TEST_ASSERT_FALSE(fc_channel_summary_seen(model, 0));
    TEST_ASSERT_TRUE(fc_channel_summary_seen(model, 1));
    TEST_ASSERT_TRUE(fc_channel_summary_seen(model, 2));
    TEST_ASSERT_FALSE(fc_channel_summary_seen(model, 3));

    fc_channel_summary_table table1 = { 0 };
    fc_channel_summary_table table2 = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_channel_summary_materialize_sparse(model, 1, &table1));
    TEST_ASSERT_EQUAL_INT(0, fc_channel_summary_materialize_sparse(model, 2, &table2));
    TEST_ASSERT_EQUAL_INT(-1, fc_channel_summary_materialize_sparse(model, 3, &table2));

    TEST_ASSERT_EQUAL_size_t(1, table1.n);
    TEST_ASSERT_EQUAL_size_t(1, table2.n);
    TEST_ASSERT_EQUAL_size_t(5, table1.rows[0].channel);
    TEST_ASSERT_EQUAL_size_t(7, table2.rows[0].channel);

    fc_channel_summary_table_destroy(&table1);
    fc_channel_summary_table_destroy(&table2);
    fc_channel_summary_destroy(model);
}

FC_TEST(channel_summary_suite, test_channel_summary_retains_unknown_mux_and_skips_invalid_values)
{
    fc_channel_summary *model = fc_channel_summary_create(1);
    TEST_ASSERT_NOT_NULL(model);

    add_channel_test_read(model, 0, 0, 1, 10, 10.0f);
    add_channel_test_read(model, 0, FC_CHANNEL_SUMMARY_MAX_CHANNEL + 1, 1, 10, 10.0f);
    add_channel_test_read(model, 0, 1, 0, 10, 10.0f);
    add_channel_test_read(model, 0, 1, FC_CHANNEL_SUMMARY_MAX_MUX + 1, 10, 10.0f);

    fc_channel_summary_skip_counts skipped = fc_channel_summary_skipped_counts(model);
    TEST_ASSERT_EQUAL_size_t(1, skipped.missing_channel);
    TEST_ASSERT_EQUAL_size_t(1, skipped.invalid_channel);
    TEST_ASSERT_EQUAL_size_t(0, skipped.missing_mux);
    TEST_ASSERT_EQUAL_size_t(1, skipped.invalid_mux);

    fc_channel_summary_table table = { 0 };
    TEST_ASSERT_EQUAL_INT(0, fc_channel_summary_materialize_sparse(model, 0, &table));
    TEST_ASSERT_EQUAL_size_t(1, table.n);
    TEST_ASSERT_EQUAL_size_t(1, table.rows[0].channel);
    TEST_ASSERT_EQUAL_size_t(0, table.rows[0].mux);
    TEST_ASSERT_EQUAL_size_t(1, table.rows[0].read_count);
    TEST_ASSERT_EQUAL_size_t(10, table.rows[0].base_count);

    fc_channel_summary_table_destroy(&table);
    fc_channel_summary_destroy(model);
}

FC_TEST(channel_summary_suite, test_channel_summary_reads_bam_mx_tag_into_meta)
{
    fc_read *read = make_channel_test_read("mx-tag");
    TEST_ASSERT_NOT_NULL(read);
    TEST_ASSERT_TRUE(bam_aux_update_int(read->bam, "ch", 42) >= 0);
    TEST_ASSERT_TRUE(bam_aux_update_int(read->bam, "mx", 3) >= 0);

    fc_read_reset_from_bam(read, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_size_t(42, read->meta->channel);
    TEST_ASSERT_EQUAL_size_t(3, read->meta->mux);

    fc_read_destroy(read);
}
