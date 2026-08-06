#ifndef FASTCAT_TEST_FS_UTILS_H
#define FASTCAT_TEST_FS_UTILS_H

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#endif
#include "core/fc_common.h"

#define FC_TEST_TEMP_PATH_CAPACITY 512

/*
 * The POSIX mkdtemp/mkstemp interfaces are not supplied by MinGW's UCRT
 * headers.  Tests use these helpers so their temporary fixtures work on both
 * Windows and POSIX hosts without changing the production code paths.
 */
static inline char *fc_test_mkdtemp(char *template_path, size_t template_size)
{
#ifdef _WIN32
    char temp_dir[MAX_PATH];
    char generated_path[MAX_PATH];
    DWORD temp_dir_len = GetTempPathA(sizeof(temp_dir), temp_dir);
    if (temp_dir_len == 0 || temp_dir_len >= sizeof(temp_dir)
        || GetTempFileNameA(temp_dir, "fct", 0, generated_path) == 0
        || !DeleteFileA(generated_path)
        || !CreateDirectoryA(generated_path, NULL)) {
        errno = EIO;
        return NULL;
    }
    size_t path_len = strlen(generated_path) + 1;
    if (path_len > template_size) {
        RemoveDirectoryA(generated_path);
        errno = ENAMETOOLONG;
        return NULL;
    }
    memcpy(template_path, generated_path, path_len);
    return template_path;
#else
    (void)template_size;
    return mkdtemp(template_path);
#endif
}

static inline int fc_test_mkstemp(char *template_path, size_t template_size)
{
#ifdef _WIN32
    char temp_dir[MAX_PATH];
    char generated_path[MAX_PATH];
    DWORD temp_dir_len = GetTempPathA(sizeof(temp_dir), temp_dir);
    if (temp_dir_len == 0 || temp_dir_len >= sizeof(temp_dir)
        || GetTempFileNameA(temp_dir, "fct", 0, generated_path) == 0) {
        errno = EIO;
        return -1;
    }
    size_t path_len = strlen(generated_path) + 1;
    if (path_len > template_size) {
        DeleteFileA(generated_path);
        errno = ENAMETOOLONG;
        return -1;
    }
    memcpy(template_path, generated_path, path_len);
    return _open(template_path, _O_RDWR | _O_BINARY);
#else
    (void)template_size;
    return mkstemp(template_path);
#endif
}

static inline int fc_test_mkdir(const char *path)
{
#ifdef _WIN32
    return mkdir(path);
#else
    return mkdir(path, 0700);
#endif
}

// Recursively removes a file tree rooted at `path`.
// Returns 0 on success (or when path does not exist), -1 on failure.
static inline int fc_test_remove_tree(const char *path) {
    return fc_remove_tree(path);
}

#endif
