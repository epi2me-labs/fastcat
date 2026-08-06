#include "unity.h"
#include "main.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "core/fc_common.h"
#include "core/fc_input_list.h"
#include "test_fs_utils.h"

static void make_file(const char *path, const char *content)
{
    FILE *fp = fopen(path, "w");
    TEST_ASSERT_NOT_NULL(fp);
    if (content) {
        fputs(content, fp);
    }
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
}

static fc_input_resolve_status resolve(
    char **paths, size_t n_paths, int recurse, bool fofn, bool filter_kind,
    fc_file_kind kind, fc_resolved_inputs *inputs)
{
    fc_path_list sources = { .paths = paths, .n_paths = n_paths };
    return fc_resolve_input_paths(&sources, recurse, fofn, filter_kind, kind, inputs);
}

FC_TEST(iter_collect_suite, test_resolve_named_fastq_and_bam_inputs)
{
    char tmpdir[] = "/tmp/fc_input_list_namedXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));
    char *fastq = fc_join_path(tmpdir, "reads.fastq");
    char *bam = fc_join_path(tmpdir, "reads.bam");
    make_file(fastq, "x");
    make_file(bam, "x");

    char *fastq_paths[] = { fastq, fastq };
    fc_resolved_inputs inputs = { 0 };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(fastq_paths, 2, 0, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    char *bam_paths[] = { bam, bam };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(bam_paths, 2, 0, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    char *mixed_paths[] = { fastq, bam };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(mixed_paths, 2, 0, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ | FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(0, fc_test_remove_tree(tmpdir));
    free(fastq);
    free(bam);
}

FC_TEST(iter_collect_suite, test_resolve_directories_collects_recognised_types)
{
    char tmpdir[] = "/tmp/fc_input_list_directoryXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));
    char *nested = fc_join_path(tmpdir, "nested");
    TEST_ASSERT_EQUAL_INT(0, mkdir(nested, 0700));
    char *top = fc_join_path(tmpdir, "top.fastq");
    char *deep = fc_join_path(nested, "deep.fq");
    char *junk = fc_join_path(tmpdir, "ignore.txt");
    make_file(top, "x");
    make_file(deep, "x");
    make_file(junk, "x");

    char *roots[] = { tmpdir };
    fc_resolved_inputs inputs = { 0 };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, 1, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(1, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, -1, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    char *bam = fc_join_path(tmpdir, "mixed.bam");
    make_file(bam, "x");
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, -1, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ | FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(3, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, -1, false, true, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ | FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(0, fc_test_remove_tree(tmpdir));
    free(top);
    free(deep);
    free(junk);
    free(bam);
    free(nested);
}

FC_TEST(iter_collect_suite, test_resolve_bam_directory_with_and_without_recursion)
{
    char tmpdir[] = "/tmp/fc_input_list_bam_directoryXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));
    char *nested = fc_join_path(tmpdir, "nested");
    TEST_ASSERT_EQUAL_INT(0, mkdir(nested, 0700));
    char *top = fc_join_path(tmpdir, "top.bam");
    char *deep = fc_join_path(nested, "deep.sam");
    make_file(top, "x");
    make_file(deep, "x");

    char *roots[] = { tmpdir };
    fc_resolved_inputs inputs = { 0 };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, 1, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(1, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(roots, 1, -1, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(0, fc_test_remove_tree(tmpdir));
    free(top);
    free(deep);
    free(nested);
}

FC_TEST(iter_collect_suite, test_resolve_fofn_infers_homogeneous_paths)
{
    char tmpdir[] = "/tmp/fc_input_list_fofnXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));
    char *first = fc_join_path(tmpdir, "first.fastq");
    char *second = fc_join_path(tmpdir, "second.fq.gz");
    char *fofn = fc_join_path(tmpdir, "inputs.fofn");
    make_file(first, "x");
    make_file(second, "x");
    FILE *fp = fopen(fofn, "w");
    TEST_ASSERT_NOT_NULL(fp);
    fprintf(fp, "%s\n%s\n", first, second);
    fclose(fp);

    char *sources[] = { fofn };
    fc_resolved_inputs inputs = { 0 };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(sources, 1, 0, true, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    char *bam = fc_join_path(tmpdir, "mixed.bam");
    make_file(bam, "x");
    fp = fopen(fofn, "w");
    TEST_ASSERT_NOT_NULL(fp);
    fprintf(fp, "%s\n%s\n", first, bam);
    fclose(fp);
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(sources, 1, 0, true, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ | FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(2, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(sources, 1, 0, true, true, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ | FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(1, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);

    TEST_ASSERT_EQUAL_INT(0, fc_test_remove_tree(tmpdir));
    free(first);
    free(second);
    free(bam);
    free(fofn);
}

FC_TEST(iter_collect_suite, test_resolve_stdin_requires_explicit_format)
{
    char *stdin_path[] = { "-" };
    fc_resolved_inputs inputs = { 0 };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_FORMAT_REQUIRED,
        resolve(stdin_path, 1, 0, false, false, FC_FILE_FASTQ, &inputs));
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_OK,
        resolve(stdin_path, 1, 0, false, true, FC_FILE_BAM, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_BAM, inputs.found_kinds);
    TEST_ASSERT_EQUAL_UINT(1, inputs.paths.n_paths);
    fc_resolved_inputs_free(&inputs);
}

FC_TEST(iter_collect_suite, test_resolve_reports_specific_rejection_statuses)
{
    char tmpdir[] = "/tmp/fc_input_list_statusXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));
    char *fastq = fc_join_path(tmpdir, "reads.fastq");
    char *unknown = fc_join_path(tmpdir, "reads.txt");
    char *empty_dir = fc_join_path(tmpdir, "empty");
    char *missing = fc_join_path(tmpdir, "missing.bam");
    make_file(fastq, "x");
    make_file(unknown, "x");
    TEST_ASSERT_EQUAL_INT(0, mkdir(empty_dir, 0700));

    fc_resolved_inputs inputs = { 0 };
    char *fastq_path[] = { fastq };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_EMPTY,
        resolve(fastq_path, 1, 0, false, true, FC_FILE_BAM, &inputs));
    TEST_ASSERT_EQUAL_UINT(FC_FILE_KINDS_FASTQ, inputs.found_kinds);
    fc_resolved_inputs_free(&inputs);

    char *unknown_path[] = { unknown };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_UNKNOWN_FORMAT,
        resolve(unknown_path, 1, 0, false, false, FC_FILE_FASTQ, &inputs));

    char *empty_path[] = { empty_dir };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_EMPTY,
        resolve(empty_path, 1, -1, false, false, FC_FILE_FASTQ, &inputs));

    char *missing_path[] = { missing };
    TEST_ASSERT_EQUAL_INT(FC_INPUT_RESOLVE_IO_ERROR,
        resolve(missing_path, 1, 0, false, false, FC_FILE_FASTQ, &inputs));

    TEST_ASSERT_EQUAL_INT(0, fc_test_remove_tree(tmpdir));
    free(fastq);
    free(unknown);
    free(empty_dir);
    free(missing);
}
