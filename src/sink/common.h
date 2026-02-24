#ifndef FASTCAT_SINK_COMMON_H
#define FASTCAT_SINK_COMMON_H

#include <stddef.h>

// shared helpers declared in src/common.h
// Returns "unclassified" for barcode 0, otherwise writes "barcode%04zu" into `buf`
// and returns `buf`. Returns NULL if `buf` is missing/too small for formatted output.
const char *fc_sink_barcode_stem(size_t barcode, char *buf, size_t buf_size);
// Returns a newly allocated path "<root>/<barcode-dir>" where barcode 0 maps to "unclassified"
// and others to "barcode%04zu". Ensures the barcode directory exists. Caller must free.
char *fc_sink_ensure_barcode_dir(const char *root, size_t barcode);

#endif
