#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "unity.h"
#include "main.h"

#include "core/fc_common.h"
#include "core/fc_read.h"
#include "core/fc_sinks.h"

static fc_read *read_sink_test_read(uint16_t flag)
{
    fc_read *read = fc_read_init();
    TEST_ASSERT_NOT_NULL(read);

    const char *seq = "ACGT";
    const char *qual = "IIII";
    TEST_ASSERT_TRUE(bam_set1(
        read->bam,
        5,
        "read",
        flag,
        (flag & BAM_FUNMAP) ? -1 : 0,
        -1,
        0,
        0,
        NULL,
        -1,
        -1,
        0,
        4,
        seq,
        qual,
        0) >= 0);
    fc_read_reset_from_bam(read, NULL, NULL, NULL);
    read->meta->runid = "secondary-run";
    read->meta->basecaller = "secondary-basecaller";
    return read;
}

static void make_temp_dir(char *path, size_t path_size)
{
    TEST_ASSERT_TRUE(path_size >= sizeof("/tmp/fastcat_read_sinks_XXXXXX"));
    FC_SNPRINTF_CHECK(path, path_size, "/tmp/fastcat_read_sinks_XXXXXX");
    TEST_ASSERT_NOT_NULL(mkdtemp(path));
}

static size_t file_line_count(const char *path)
{
    FILE *fp = fopen(path, "r");
    TEST_ASSERT_NOT_NULL(fp);
    size_t lines = 0;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        ++lines;
    }
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
    return lines;
}

static bool file_contains(const char *path, const char *needle)
{
    FILE *fp = fopen(path, "r");
    TEST_ASSERT_NOT_NULL(fp);
    char buffer[4096] = { 0 };
    size_t n = fread(buffer, 1, sizeof(buffer) - 1, fp);
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
    buffer[n] = '\0';
    return strstr(buffer, needle) != NULL;
}

FC_TEST(read_sinks_suite, test_output_and_perread_sinks_ignore_additional_alignments)
{
    char output_dir[64];
    make_temp_dir(output_dir, sizeof(output_dir));
    fc_read *read = read_sink_test_read(BAM_FSECONDARY | BAM_FUNMAP);
    fc_read_ctx ctx = { 0 };

    fc_sink output = fc_output_sink_create(&(fc_output_sink_opts) {
        .output_root = output_dir,
        .barcode_count = 1
    });
    TEST_ASSERT_NOT_NULL(output.ctx);
    output.add_read(output.ctx, read, &ctx, "input.bam", NULL);
    output.destroy(output.ctx);

    char path[128];
    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/unclassified", output_dir);
    TEST_ASSERT_EQUAL_INT(-1, access(path, F_OK));

    fc_sink perread = fc_perread_sink_create(&(fc_perread_sink_opts) {
        .output_root = output_dir,
        .barcode_count = 1
    });
    TEST_ASSERT_NOT_NULL(perread.ctx);
    perread.add_read(perread.ctx, read, &ctx, "input.bam", NULL);
    perread.destroy(perread.ctx);

    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/per_read.tsv", output_dir);
    TEST_ASSERT_EQUAL_INT(-1, access(path, F_OK));

    fc_read_destroy(read);
    TEST_ASSERT_EQUAL_INT(0, rmdir(output_dir));
}

FC_TEST(read_sinks_suite, test_metadata_counters_ignore_additional_alignments)
{
    char output_dir[64];
    make_temp_dir(output_dir, sizeof(output_dir));
    fc_read *read = read_sink_test_read(BAM_FSUPPLEMENTARY | BAM_FUNMAP);
    fc_read_ctx ctx = { 0 };

    fc_sink counter = fc_counter_sink_create(&(fc_counter_sink_opts) {
        .output_root = output_dir,
        .barcode_count = 1
    });
    TEST_ASSERT_NOT_NULL(counter.ctx);
    counter.begin_file(counter.ctx, "input.bam");
    counter.add_read(counter.ctx, read, &ctx, "input.bam", NULL);
    counter.end_file(counter.ctx, "input.bam", &(fc_file_stats) { 0 });
    counter.destroy(counter.ctx);

    char path[128];
    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/runids.tsv", output_dir);
    TEST_ASSERT_EQUAL_size_t(1, file_line_count(path));
    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/basecallers.tsv", output_dir);
    TEST_ASSERT_EQUAL_size_t(1, file_line_count(path));

    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/runids.tsv", output_dir);
    TEST_ASSERT_EQUAL_INT(0, remove(path));
    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/basecallers.tsv", output_dir);
    TEST_ASSERT_EQUAL_INT(0, remove(path));
    fc_read_destroy(read);
    TEST_ASSERT_EQUAL_INT(0, rmdir(output_dir));
}

FC_TEST(read_sinks_suite, test_summary_json_sink_ignores_additional_alignments)
{
    char output_dir[64];
    make_temp_dir(output_dir, sizeof(output_dir));
    fc_read *read = read_sink_test_read(BAM_FSECONDARY | BAM_FUNMAP);
    fc_read_ctx ctx = { 0 };
    fc_stats_engine *engine = fc_stats_engine_create(false, false, 0.0f, 0.0f, false);
    fc_stats_engine *engines[] = { engine };
    fc_run_summary *run_summary = fc_run_summary_create(1);
    fc_channel_summary *channel_summary = fc_channel_summary_create(1);
    fc_length_qscore_2d *length_qscore = fc_length_qscore_2d_create(1, 64, 33, 1);
    fc_speed_qscore_2d *speed_qscore = fc_speed_qscore_2d_create(1, 64, 33, 1);
    fc_time_bin *time_bin = fc_time_bin_create(1, 10.0, 64);
    TEST_ASSERT_NOT_NULL(engine);
    TEST_ASSERT_NOT_NULL(run_summary);
    TEST_ASSERT_NOT_NULL(channel_summary);
    TEST_ASSERT_NOT_NULL(length_qscore);
    TEST_ASSERT_NOT_NULL(speed_qscore);
    TEST_ASSERT_NOT_NULL(time_bin);

    fc_sink summary_json = fc_summary_json_sink_create(&(fc_summary_json_sink_opts) {
        .output_root = output_dir,
        .engines = engines,
        .engine_count = 1,
        .run_summary = run_summary,
        .channel_summary = channel_summary,
        .length_qscore = length_qscore,
        .speed_qscore = speed_qscore,
        .time_bin = time_bin
    });
    TEST_ASSERT_NOT_NULL(summary_json.ctx);
    summary_json.begin_file(summary_json.ctx, "input.bam");
    summary_json.add_read(summary_json.ctx, read, &ctx, "input.bam", NULL);
    summary_json.end_file(summary_json.ctx, "input.bam", &(fc_file_stats) { 0 });
    summary_json.finish(summary_json.ctx);
    summary_json.destroy(summary_json.ctx);

    char path[128];
    FC_SNPRINTF_CHECK(path, sizeof(path), "%s/summary.json", output_dir);
    TEST_ASSERT_FALSE(file_contains(path, "secondary-run"));
    TEST_ASSERT_EQUAL_INT(0, remove(path));

    fc_time_bin_destroy(time_bin);
    fc_speed_qscore_2d_destroy(speed_qscore);
    fc_length_qscore_2d_destroy(length_qscore);
    fc_channel_summary_destroy(channel_summary);
    fc_run_summary_destroy(run_summary);
    fc_stats_engine_destroy(engine);
    fc_read_destroy(read);
    TEST_ASSERT_EQUAL_INT(0, rmdir(output_dir));
}

FC_TEST(read_sinks_suite, test_quality_heatmap_models_ignore_additional_alignments)
{
    fc_read *read = read_sink_test_read(BAM_FSECONDARY | BAM_FUNMAP);
    read->meta->mean_q = 20.0f;
    read->meta->duration_s = 2.0;
    fc_length_qscore_2d *length_qscore = fc_length_qscore_2d_create(1, 64, 33, 1);
    fc_speed_qscore_2d *speed_qscore = fc_speed_qscore_2d_create(1, 64, 33, 1);
    fc_accuracy_qscore_2d *accuracy_qscore = fc_accuracy_qscore_2d_create(1, 64, 33, 1);
    TEST_ASSERT_NOT_NULL(length_qscore);
    TEST_ASSERT_NOT_NULL(speed_qscore);
    TEST_ASSERT_NOT_NULL(accuracy_qscore);

    fc_length_qscore_2d_add_read(length_qscore, 0, read);
    fc_speed_qscore_2d_add_read(speed_qscore, 0, read);
    fc_accuracy_qscore_2d_add_read(accuracy_qscore, 0, read, 99.0f);
    TEST_ASSERT_FALSE(fc_length_qscore_2d_seen(length_qscore, 0));
    TEST_ASSERT_FALSE(fc_speed_qscore_2d_seen(speed_qscore, 0));
    TEST_ASSERT_FALSE(fc_accuracy_qscore_2d_seen(accuracy_qscore, 0));

    read->bam->core.flag = BAM_FUNMAP;
    fc_length_qscore_2d_add_read(length_qscore, 0, read);
    fc_speed_qscore_2d_add_read(speed_qscore, 0, read);
    fc_accuracy_qscore_2d_add_read(accuracy_qscore, 0, read, 99.0f);
    TEST_ASSERT_TRUE(fc_length_qscore_2d_seen(length_qscore, 0));
    TEST_ASSERT_TRUE(fc_speed_qscore_2d_seen(speed_qscore, 0));
    TEST_ASSERT_TRUE(fc_accuracy_qscore_2d_seen(accuracy_qscore, 0));

    fc_accuracy_qscore_2d_destroy(accuracy_qscore);
    fc_speed_qscore_2d_destroy(speed_qscore);
    fc_length_qscore_2d_destroy(length_qscore);
    fc_read_destroy(read);
}
