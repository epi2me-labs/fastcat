#ifndef FASTCAT_TEST_FS_UTILS_H
#define FASTCAT_TEST_FS_UTILS_H

#include "core/fc_common.h"

// Recursively removes a file tree rooted at `path`.
// Returns 0 on success (or when path does not exist), -1 on failure.
static int fc_test_remove_tree(const char *path) {
    return fc_remove_tree(path);
}

#endif
