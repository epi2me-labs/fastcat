#ifndef FASTCAT_INPUT_LIST_H
#define FASTCAT_INPUT_LIST_H

#include <stddef.h>

#include "fc_iter.h"

// Read newline-delimited file paths from stdin in --lofn mode.
// Validates each path as an existing file of expected_kind (no directories,
// no "-" entries). Returns 0 on success and fills out_paths.
// Caller must free out_paths with fc_path_list_free.
int fc_read_paths_from_stdin(fc_file_kind expected_kind, fc_path_list *out_paths);

#endif
