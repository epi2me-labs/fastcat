#ifndef FASTCAT_INPUT_LIST_H
#define FASTCAT_INPUT_LIST_H

#include <stddef.h>

// Read newline-delimited paths from stdin. Returns a heap array of strings and
// sets out_n to the count. Caller must free with fc_free_path_list.
char **fc_read_paths_from_stdin(size_t *out_n);

// Free a path list allocated by fc_read_paths_from_stdin.
void fc_free_path_list(char **paths, size_t n);

#endif
