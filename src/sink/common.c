#include "sink/common.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

char *fc_sink_strip_path(const char *input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (len > 0 && input[len - 1] == '/') len--;
    char *output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
}

char *fc_sink_join_path(const char *dir, const char *leaf) {
    if (!dir || !leaf) return NULL;
    size_t len = strlen(dir);
    size_t llen = strlen(leaf);
    bool need_sep = len > 0 && dir[len - 1] != '/';
    size_t total = len + need_sep + llen + 1;
    char *out = calloc(total, 1);
    if (!out) return NULL;
    strcpy(out, dir);
    if (need_sep) strcat(out, "/");
    strcat(out, leaf);
    return out;
}
