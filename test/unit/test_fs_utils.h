#ifndef FASTCAT_TEST_FS_UTILS_H
#define FASTCAT_TEST_FS_UTILS_H

#include <dirent.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "core/fc_common.h"

// Recursively removes a file tree rooted at `path`.
// Returns 0 on success (or when path does not exist), -1 on failure.
static int fc_test_remove_tree(const char *path) {
    if (!path) {
        errno = EINVAL;
        return -1;
    }

    struct stat st = {0};
    if (lstat(path, &st) != 0) {
        if (errno == ENOENT) return 0;
        return -1;
    }

    if (!S_ISDIR(st.st_mode)) {
        return unlink(path);
    }

    DIR *dir = opendir(path);
    if (!dir) return -1;

    int status = 0;
    int saved_errno = 0;
    struct dirent *ent = NULL;
    while ((ent = readdir(dir)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
            continue;
        }

        size_t clen = strlen(path) + 1 + strlen(ent->d_name) + 1;
        char *child = malloc(clen);
        if (!child) {
            status = -1;
            saved_errno = ENOMEM;
            break;
        }
        FC_SNPRINTF_CHECK(child, clen, "%s/%s", path, ent->d_name);
        if (fc_test_remove_tree(child) != 0 && status == 0) {
            status = -1;
            saved_errno = errno;
        }
        free(child);
    }

    if (closedir(dir) != 0 && status == 0) {
        status = -1;
        saved_errno = errno;
    }

    if (rmdir(path) != 0 && status == 0) {
        status = -1;
        saved_errno = errno;
    }

    if (status != 0 && saved_errno != 0) {
        errno = saved_errno;
    }
    return status;
}

#endif
