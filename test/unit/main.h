#ifndef TEST_UNIT_MAIN_H
#define TEST_UNIT_MAIN_H

#include "core/fc_bam_rg_cache.h"

typedef void (*test_fn)(void);

void fc_register_test(
    const char *suite_name,
    const char *test_name,
    const char *file,
    int line,
    test_fn fn);
void fc_run_registered_tests(void);
fc_bam_rg_cache *fc_test_rg_cache(void);

// Each test uses this macro instead of a plain `static void test_foo(void)`.
// It defines the real test function and a constructor that registers the
// function pointer together with its original source file and line number.
// `test/unit/main.c` then sorts the collected tests by suite name and test
// name and dispatches them via Unity directly, restoring sensible file/line
// reporting without needing per-file RUN_TEST(...) lists.
#define FC_TEST(suite_name, test_name) \
    static void test_name(void); \
    static void fc_register_test_##suite_name##_##test_name(void) __attribute__((constructor)); \
    static void fc_register_test_##suite_name##_##test_name(void) \
    { \
        fc_register_test(#suite_name, #test_name, __FILE__, __LINE__, test_name); \
    } \
    static void test_name(void)

#endif
