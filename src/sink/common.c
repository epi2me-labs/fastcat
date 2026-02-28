#include "sink/common.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common.h" // ensure_parent_dir_exists, path helpers

const char *fc_sink_barcode_stem(size_t barcode, char *buf, size_t buf_size)
{
    if (barcode == 0) {
        return "unclassified";
    }
    if (!buf || buf_size == 0) {
        return NULL;
    }
    FC_SNPRINTF_CHECK(buf, buf_size, "barcode%04zu", barcode);
    return buf;
}

char *fc_sink_ensure_barcode_dir(const char *root, size_t barcode)
{
    if (!root) {
        return NULL;
    }
    char barc[32] = { 0 };
    const char *sub = fc_sink_barcode_stem(barcode, barc, sizeof(barc));
    if (!sub) {
        return NULL;
    }
    char *dir = fc_join_path(root, sub);
    // mkdir_p ensures the barcode directory itself exists (not just its parent)
    if (mkdir_p(dir) != 0) {
        free(dir);
        return NULL;
    }
    return dir;
}
