#ifndef FASTCAT_INPUT_LIST_H
#define FASTCAT_INPUT_LIST_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    FC_FILE_FASTQ = 0,
    FC_FILE_BAM
} fc_file_kind;

enum {
    FC_FILE_KINDS_NONE = 0,
    FC_FILE_KINDS_FASTQ = 1 << FC_FILE_FASTQ,
    FC_FILE_KINDS_BAM = 1 << FC_FILE_BAM
};

typedef struct {
    char **paths;
    size_t n_paths;
} fc_path_list;

typedef struct {
    fc_path_list paths;
    fc_file_kind kind;
    unsigned found_kinds;
} fc_resolved_inputs;

typedef enum {
    FC_INPUT_RESOLVE_OK = 0,
    FC_INPUT_RESOLVE_INVALID_ARGUMENT,
    FC_INPUT_RESOLVE_IO_ERROR,
    FC_INPUT_RESOLVE_EMPTY,
    FC_INPUT_RESOLVE_FORMAT_REQUIRED,
    FC_INPUT_RESOLVE_UNKNOWN_FORMAT,
    FC_INPUT_RESOLVE_MIXED_FORMAT,
    FC_INPUT_RESOLVE_FOFN_ERROR
} fc_input_resolve_status;

bool fc_path_matches_kind(const char *path, fc_file_kind kind);
void fc_path_list_free(fc_path_list *paths);
void fc_resolved_inputs_free(fc_resolved_inputs *inputs);

// Expand positional paths, directories, or FOFN sources. Recognised FASTQ and
// BAM types are recorded in found_kinds. When filter_kind is true, only files
// of requested_kind are returned. stdin data requires a filter kind.
fc_input_resolve_status fc_resolve_input_paths(
    const fc_path_list *sources,
    int recurse,
    bool stdin_is_path_list,
    bool filter_kind,
    fc_file_kind requested_kind,
    fc_resolved_inputs *out_inputs);
const char *fc_input_resolve_status_string(fc_input_resolve_status status);

#endif
