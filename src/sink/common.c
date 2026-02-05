#include "sink/common.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

#include "../common.h"  // ensure_parent_dir_exists, path helpers

char *fc_sink_ensure_barcode_dir(const char *root, size_t barcode) {
    if (!root) return NULL;
    const char *sub = (barcode == 0) ? "unclassified" : NULL;
    char barc[32];
    if (!sub) {
        snprintf(barc, sizeof(barc), "barcode%04zu", barcode);
        sub = barc;
    }
    char *dir = fc_join_path(root, sub);
    // mkdir_p ensures the barcode directory itself exists (not just its parent)
    if (mkdir_p(dir) != 0) {
        free(dir);
        return NULL;
    }
    return dir;
}
