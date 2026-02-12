#include "fc_input_list.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fc_log.h"

char **fc_read_paths_from_stdin(size_t *out_n) {
    if (out_n) *out_n = 0;
    char **paths = NULL;
    char *line = NULL;
    size_t cap = 0;
    ssize_t nchr = 0;
    while ((nchr = getline(&line, &cap, stdin)) != -1) {
        line[strcspn(line, "\r\n")] = 0;
        if (line[0] == '\0') continue;
        char *copy = strdup(line);
        if (!copy) {
            FC_LOG_ERROR("Out of memory reading input list.");
            break;
        }
        char **tmp = realloc(paths, (*out_n + 1) * sizeof(char *));
        if (!tmp) {
            FC_LOG_ERROR("Out of memory growing input list.");
            free(copy);
            break;
        }
        paths = tmp;
        paths[*out_n] = copy;
        (*out_n)++;
    }
    free(line);
    FC_LOG_INFO("Read %zu paths from stdin.", *out_n);
    return paths;
}

void fc_free_path_list(char **paths, size_t n) {
    if (!paths) return;
    for (size_t i = 0; i < n; ++i) {
        free(paths[i]);
    }
    free(paths);
}
