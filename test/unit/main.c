#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "core/fc_log.h"
#include "main.h"

typedef struct {
    const char *suite_name;
    const char *test_name;
    const char *file;
    int line;
    test_fn fn;
} test_entry;

static test_entry *tests = NULL;
static size_t test_count = 0;
static size_t test_capacity = 0;

void setUp(void) {}
void tearDown(void) {}

void fc_register_test(
    const char *suite_name,
    const char *test_name,
    const char *file,
    int line,
    test_fn fn)
{
    if (NULL == suite_name || NULL == test_name || NULL == file || line <= 0 || NULL == fn) {
        return;
    }
    if (test_count == test_capacity) {
        size_t next_capacity = (0 == test_capacity) ? 32 : test_capacity * 2;
        test_entry *next = realloc(tests, next_capacity * sizeof(*next));
        if (NULL == next) {
            fprintf(stderr, "failed to register test %s/%s\n", suite_name, test_name);
            abort();
        }
        tests = next;
        test_capacity = next_capacity;
    }
    tests[test_count].suite_name = suite_name;
    tests[test_count].test_name = test_name;
    tests[test_count].file = file;
    tests[test_count].line = line;
    tests[test_count].fn = fn;
    test_count++;
}

static int compare_test_entries(const void *lhs, const void *rhs)
{
    const test_entry *left = lhs;
    const test_entry *right = rhs;
    int suite_cmp = strcmp(left->suite_name, right->suite_name);
    if (suite_cmp != 0) {
        return suite_cmp;
    }
    return strcmp(left->test_name, right->test_name);
}

void fc_run_registered_tests(void)
{
    qsort(tests, test_count, sizeof(*tests), compare_test_entries);
    const char *current_suite = NULL;
    for (size_t i = 0; i < test_count; i++) {
        if (NULL == current_suite || 0 != strcmp(current_suite, tests[i].suite_name)) {
            current_suite = tests[i].suite_name;
            fprintf(stderr, "\nRunning test suite: %s\n", current_suite);
        }
        UnitySetTestFile(tests[i].file);
        UnityDefaultTestRun(tests[i].fn, tests[i].test_name, tests[i].line);
    }
}

int main(void)
{
    fc_log_config log_cfg = {
        .level = FC_LOG_DEBUG,
        .module = "UNIT_TEST",
        .target = "/dev/null",
        .color = FC_LOG_COLOR_OFF
    };
    fc_log_init(&log_cfg);

    UNITY_BEGIN();
    fc_run_registered_tests();
    int rc = UNITY_END();

    fc_log_destroy();
    free(tests);
    return rc;
}
