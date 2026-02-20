#include "fc_input_list.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "fc_log.h"

int fc_read_paths_from_stdin(fc_file_kind expected_kind, fc_path_list *out_paths) {
    if (!out_paths) return EINVAL;
    out_paths->paths = NULL;
    out_paths->n_paths = 0;

    int status = 0;
    char **paths = NULL;
    char *line = NULL;
    size_t cap = 0;
    ssize_t nchr = 0;
    while ((nchr = getline(&line, &cap, stdin)) != -1) {
        struct stat st = {0};
        line[strcspn(line, "\r\n")] = 0;
        if (line[0] == '\0') continue;
        if (strcmp(line, "-") == 0) {
            FC_LOG_ERROR("'-' is not valid in --lofn input; provide file paths only.");
            status = (status == 0) ? EINVAL : status;
            continue;
        }
        if (stat(line, &st) != 0) {
            FC_LOG_ERROR("could not stat '%s': %s", line, strerror(errno));
            status = (status == 0) ? errno : status;
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            FC_LOG_ERROR("path list input '%s' is a directory; expected a file path.", line);
            status = (status == 0) ? EISDIR : status;
            continue;
        }
        if (!fc_path_matches_kind(line, expected_kind)) {
            FC_LOG_ERROR("Input '%s' does not match required input type for this command.", line);
            status = (status == 0) ? EINVAL : status;
            continue;
        }
        char *copy = strdup(line);
        if (!copy) {
            FC_LOG_ERROR("Out of memory reading input list.");
            status = ENOMEM;
            break;
        }
        char **tmp = realloc(paths, (out_paths->n_paths + 1) * sizeof(char *));
        if (!tmp) {
            FC_LOG_ERROR("Out of memory growing input list.");
            free(copy);
            status = ENOMEM;
            break;
        }
        paths = tmp;
        paths[out_paths->n_paths] = copy;
        out_paths->n_paths++;
    }
    free(line);
    FC_LOG_INFO("Read %zu paths from stdin.", out_paths->n_paths);
    if (status == 0 && out_paths->n_paths == 0) {
        status = EINVAL;
    }
    if (status != 0) {
        fc_path_list tmp = {.paths = paths, .n_paths = out_paths->n_paths};
        fc_path_list_free(&tmp);
        out_paths->n_paths = 0;
        return status;
    }
    out_paths->paths = paths;
    return 0;
}
