#ifndef FASTCAT_SINK_COMMON_H
#define FASTCAT_SINK_COMMON_H

#include <stddef.h>

// shared helpers declared in src/common.h
// Returns a newly allocated path "<root>/<barcode-dir>" where barcode 0 maps to "unclassified"
// and others to "barcode%04zu". Ensures the barcode directory exists. Caller must free.
char *fc_sink_ensure_barcode_dir(const char *root, size_t barcode);

#endif
