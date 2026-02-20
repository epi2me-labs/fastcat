#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#include "fc_iter.h"
#include "common.h"  // fc_join_path

static void make_file(const char *path, const char *content) {
    FILE *fp = fopen(path, "w");
    TEST_ASSERT_NOT_NULL(fp);
    if (content) fputs(content, fp);
    TEST_ASSERT_EQUAL_INT(0, fclose(fp));
}

static void test_collect_paths_filters_and_recurses_bfs(void) {
    char tmpdir[] = "/tmp/fc_iter_collectXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));

    char *sub1 = fc_join_path(tmpdir, "sub1");
    TEST_ASSERT_NOT_NULL(sub1);
    char *sub2 = fc_join_path(tmpdir, "sub2");
    TEST_ASSERT_NOT_NULL(sub2);
    char *junkdir = fc_join_path(tmpdir, "junk");
    TEST_ASSERT_NOT_NULL(junkdir);
    TEST_ASSERT_EQUAL_INT(0, mkdir(sub1, 0700));
    TEST_ASSERT_EQUAL_INT(0, mkdir(sub2, 0700));
    TEST_ASSERT_EQUAL_INT(0, mkdir(junkdir, 0700));

    // valid FASTQ files (multiple per level)
    char *f1 = fc_join_path(tmpdir, "a.fastq");
    char *f1b = fc_join_path(tmpdir, "a2.fq");
    char *f2 = fc_join_path(sub1, "b.fastq");
    char *f2b = fc_join_path(sub1, "b2.fq.gz");
    char *f3 = fc_join_path(sub2, "c.fastq.gz");
    make_file(f1, "x");
    make_file(f1b, "x");
    make_file(f2, "x");
    make_file(f2b, "x");
    make_file(f3, "x");

    // junk files and nested junk dir
    char *junkfile = fc_join_path(tmpdir, "ignore.txt");
    make_file(junkfile, "x");
    char *junkfile_sub = fc_join_path(sub1, "skip.me");
    make_file(junkfile_sub, "x"); // should be ignored due to extension
    char *nestedjunkdir = fc_join_path(junkdir, "nest");
    TEST_ASSERT_EQUAL_INT(0, mkdir(nestedjunkdir, 0700));
    char *nestedjunk = fc_join_path(nestedjunkdir, "skip.fq.bak");
    make_file(nestedjunk, "x");

    char *roots[2] = { tmpdir, NULL };
    fc_path_list input_paths = {.paths = roots, .n_paths = 1};
    fc_path_list files = {0};
    TEST_ASSERT_EQUAL_INT(0, fc_iter_collect_paths(&input_paths, 2, FC_FILE_FASTQ, false, &files));

    TEST_ASSERT_EQUAL_UINT(5, files.n_paths);
    // order should be BFS by depth
    // collect the basenames to assert membership
    int seen_fastq = 0, seen_fastq2 = 0, seen_sub_fastq = 0, seen_sub_fq_gz = 0, seen_fastq_gz = 0;
    for (size_t i = 0; i < files.n_paths; ++i) {
        if (strstr(files.paths[i], "a.fastq")) seen_fastq++;
        if (strstr(files.paths[i], "a2.fq")) seen_fastq2++;
        if (strstr(files.paths[i], "b.fastq")) seen_sub_fastq++;
        if (strstr(files.paths[i], "b2.fq.gz")) seen_sub_fq_gz++;
        if (strstr(files.paths[i], "c.fastq.gz")) seen_fastq_gz++;
    }
    TEST_ASSERT_EQUAL_INT(1, seen_fastq);
    TEST_ASSERT_EQUAL_INT(1, seen_fastq2);
    TEST_ASSERT_EQUAL_INT(1, seen_sub_fastq);
    TEST_ASSERT_EQUAL_INT(1, seen_sub_fq_gz);
    TEST_ASSERT_EQUAL_INT(1, seen_fastq_gz);

    fc_path_list_free(&files);

    // cleanup files and dirs like a good citizen
    unlink(f1); unlink(f1b); unlink(f2); unlink(f2b); unlink(f3);
    unlink(junkfile); unlink(junkfile_sub); unlink(nestedjunk);
    rmdir(nestedjunkdir);
    rmdir(junkdir);
    rmdir(sub1);
    rmdir(sub2);
    rmdir(tmpdir);
    free(f1); free(f1b); free(f2); free(f2b); free(f3);
    free(junkfile); free(junkfile_sub); free(nestedjunk); free(nestedjunkdir);
    free(sub1); free(sub2); free(junkdir);
}

static void test_collect_paths_rejects_mixed_file_types(void) {
    char tmpdir[] = "/tmp/fc_iter_collect_mixedXXXXXX";
    TEST_ASSERT_NOT_NULL(mkdtemp(tmpdir));

    char *f1 = fc_join_path(tmpdir, "a.fastq");
    char *f2 = fc_join_path(tmpdir, "b.sam");
    make_file(f1, "x");
    make_file(f2, "x");

    char *roots[2] = { tmpdir, NULL };
    fc_path_list input_paths = {.paths = roots, .n_paths = 1};
    fc_path_list files = {0};
    int rc = fc_iter_collect_paths(&input_paths, 1, FC_FILE_FASTQ, false, &files);
    TEST_ASSERT_EQUAL_INT(EINVAL, rc);
    TEST_ASSERT_NULL(files.paths);
    TEST_ASSERT_EQUAL_UINT(0, files.n_paths);

    unlink(f1);
    unlink(f2);
    rmdir(tmpdir);
    free(f1);
    free(f2);
}

void iter_collect_suite(void) {
    RUN_TEST(test_collect_paths_filters_and_recurses_bfs);
    RUN_TEST(test_collect_paths_rejects_mixed_file_types);
}
