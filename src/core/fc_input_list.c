#include "core/fc_input_list.h"

#include <dirent.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "core/fc_common.h"
#include "core/fc_log.h"

static fc_input_resolve_status read_fofn_sources(
    const fc_path_list *sources,
    fc_path_list *out_paths);

static bool has_suffix(const char *path, const char *suffix)
{
    size_t path_len = strlen(path);
    size_t suffix_len = strlen(suffix);
    return path_len >= suffix_len && strcmp(path + path_len - suffix_len, suffix) == 0;
}

bool fc_path_matches_kind(const char *path, fc_file_kind kind)
{
    if (!path) {
        return false;
    }
    if (kind == FC_FILE_FASTQ) {
        return has_suffix(path, ".fastq") || has_suffix(path, ".fastq.gz")
            || has_suffix(path, ".fq") || has_suffix(path, ".fq.gz");
    }
    return kind == FC_FILE_BAM && (has_suffix(path, ".bam") || has_suffix(path, ".sam"));
}

static int append_path(fc_path_list *list, const char *path)
{
    char *copy = fc_strdup(path, "input list path");
    list->paths = fc_realloc(list->paths, (list->n_paths + 1) * sizeof(*list->paths), "input list paths array");
    list->paths[list->n_paths] = copy;
    list->n_paths++;
    return 0;
}

void fc_path_list_free(fc_path_list *paths)
{
    if (!paths || !paths->paths) {
        return;
    }
    for (size_t i = 0; i < paths->n_paths; ++i) {
        free(paths->paths[i]);
    }
    free(paths->paths);
    paths->paths = NULL;
    paths->n_paths = 0;
}

void fc_resolved_inputs_free(fc_resolved_inputs *inputs)
{
    if (inputs) {
        fc_path_list_free(&inputs->paths);
    }
}

static bool detect_kind(const char *path, fc_file_kind *out_kind)
{
    if (fc_path_matches_kind(path, FC_FILE_FASTQ)) {
        *out_kind = FC_FILE_FASTQ;
    } else if (fc_path_matches_kind(path, FC_FILE_BAM)) {
        *out_kind = FC_FILE_BAM;
    } else {
        return false;
    }
    return true;
}

typedef struct {
    char *path;
    int depth;
    bool from_directory;
} pending_path;

static fc_input_resolve_status resolve_roots(
    const fc_path_list *sources,
    int recurse,
    bool filter_kind,
    fc_file_kind requested_kind,
    fc_resolved_inputs *out)
{
    pending_path *queue = NULL;
    size_t n_queue = 0;
    size_t queue_cap = 0;
    size_t next = 0;
    fc_input_resolve_status status = FC_INPUT_RESOLVE_OK;

    for (size_t i = 0; i < sources->n_paths; ++i) {
        if (n_queue == queue_cap) {
            size_t new_cap = queue_cap ? queue_cap * 2 : 8;
            queue = fc_realloc(queue, new_cap * sizeof(*queue), "input path queue");
            queue_cap = new_cap;
        }
        queue[n_queue++] = (pending_path) {
            .path = fc_strdup(sources->paths[i], "queued input path"),
            .depth = recurse,
            .from_directory = false
        };
    }

    while (next < n_queue && status == FC_INPUT_RESOLVE_OK) {
        pending_path current = queue[next++];
        if (strcmp(current.path, "-") == 0) {
            if (!filter_kind || sources->n_paths != 1) {
                FC_LOG_ERROR("stdin input requires --input-fmt and cannot be combined with file paths.");
                status = FC_INPUT_RESOLVE_FORMAT_REQUIRED;
            } else {
                append_path(&out->paths, current.path);
                out->found_kinds |= 1U << requested_kind;
            }
            free(current.path);
            continue;
        }

        fc_file_stat st = { 0 };
        if (fc_stat_path(current.path, &st) != 0) {
            FC_LOG_ERROR("could not stat '%s': %s", current.path, strerror(errno));
            status = FC_INPUT_RESOLVE_IO_ERROR;
            free(current.path);
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            if (current.depth != 0) {
                DIR *dir = opendir(current.path);
                if (!dir) {
                    FC_LOG_ERROR("could not open directory '%s': %s", current.path, strerror(errno));
                    status = FC_INPUT_RESOLVE_IO_ERROR;
                } else {
                    struct dirent *entry;
                    while (status == FC_INPUT_RESOLVE_OK && (entry = readdir(dir)) != NULL) {
                        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
                            continue;
                        }
                        if (n_queue == queue_cap) {
                            size_t new_cap = queue_cap ? queue_cap * 2 : 8;
                            queue = fc_realloc(queue, new_cap * sizeof(*queue), "recursive input path queue");
                            queue_cap = new_cap;
                        }
                        queue[n_queue++] = (pending_path) {
                            .path = fc_asprintf("recursive input child", "%s/%s", current.path, entry->d_name),
                            .depth = current.depth < 0 ? -1 : current.depth - 1,
                            .from_directory = true
                        };
                    }
                    closedir(dir);
                }
            }
            free(current.path);
            continue;
        }

        fc_file_kind found_kind;
        bool recognised = detect_kind(current.path, &found_kind);
        if (!recognised && current.from_directory) {
            // Preserve directory discovery semantics: unrelated files are ignored.
            free(current.path);
            continue;
        }
        if (!recognised) {
            FC_LOG_ERROR("Could not infer input format from '%s'.", current.path);
            status = FC_INPUT_RESOLVE_UNKNOWN_FORMAT;
        } else {
            out->found_kinds |= 1U << found_kind;
            if (!filter_kind || found_kind == requested_kind) {
                append_path(&out->paths, current.path);
            }
        }
        free(current.path);
    }
    for (; next < n_queue; ++next) {
        free(queue[next].path);
    }
    free(queue);
    if (status != FC_INPUT_RESOLVE_OK) {
        return status;
    }
    out->kind = requested_kind;
    return 0;
}

fc_input_resolve_status fc_resolve_input_paths(
    const fc_path_list *sources,
    int recurse,
    bool stdin_is_path_list,
    bool filter_kind,
    fc_file_kind requested_kind,
    fc_resolved_inputs *out_inputs)
{
    if (!sources || !sources->paths || sources->n_paths == 0 || !out_inputs) {
        return FC_INPUT_RESOLVE_INVALID_ARGUMENT;
    }
    memset(out_inputs, 0, sizeof(*out_inputs));
    out_inputs->kind = requested_kind;

    if (stdin_is_path_list) {
        fc_path_list listed = { 0 };
        fc_input_resolve_status read_status = read_fofn_sources(
            sources, &listed);
        if (read_status != FC_INPUT_RESOLVE_OK) {
            return read_status;
        }
        fc_input_resolve_status status = resolve_roots(
            &listed, 0, filter_kind, requested_kind, out_inputs);
        fc_path_list_free(&listed);
        if (status == FC_INPUT_RESOLVE_OK && out_inputs->paths.n_paths == 0) {
            status = FC_INPUT_RESOLVE_EMPTY;
        }
        return status;
    }

    fc_input_resolve_status status = resolve_roots(
        sources, recurse, filter_kind, requested_kind, out_inputs);
    if (status == FC_INPUT_RESOLVE_OK && out_inputs->paths.n_paths == 0) {
        status = FC_INPUT_RESOLVE_EMPTY;
    }
    if (status != FC_INPUT_RESOLVE_OK && status != FC_INPUT_RESOLVE_EMPTY) {
        fc_resolved_inputs_free(out_inputs);
    }
    return status;
}

const char *fc_input_resolve_status_string(fc_input_resolve_status status)
{
    switch (status) {
    case FC_INPUT_RESOLVE_OK:
        return "success";
    case FC_INPUT_RESOLVE_INVALID_ARGUMENT:
        return "invalid input resolver arguments";
    case FC_INPUT_RESOLVE_IO_ERROR:
        return "could not access an input path";
    case FC_INPUT_RESOLVE_EMPTY:
        return "no usable input files found";
    case FC_INPUT_RESOLVE_FORMAT_REQUIRED:
        return "standard input requires --input-fmt";
    case FC_INPUT_RESOLVE_UNKNOWN_FORMAT:
        return "could not infer an input format";
    case FC_INPUT_RESOLVE_MIXED_FORMAT:
        return "input paths have mixed formats";
    case FC_INPUT_RESOLVE_FOFN_ERROR:
        return "could not read a file-of-filenames input";
    }
    return "unknown input resolution error";
}

static fc_input_resolve_status read_paths_from_stream(
    FILE *fp,
    const char *source_name,
    fc_path_list *list)
{
    fc_input_resolve_status status = FC_INPUT_RESOLVE_OK;
    char *line = NULL;
    size_t cap = 0;
    ptrdiff_t nchr = 0;

    while ((nchr = fc_getline(&line, &cap, fp)) != -1) {
        fc_file_stat st = { 0 };
        // Remove trailing newline
        if (nchr > 0 && (line[nchr - 1] == '\n' || line[nchr - 1] == '\r')) {
            line[nchr - 1] = '\0';
            if (nchr > 1 && line[nchr - 2] == '\r') {
                line[nchr - 2] = '\0';
            }
        }
        if (line[0] == '\0') {
            continue;
        }

        if (strcmp(line, "-") == 0) {
            FC_LOG_ERROR("'-' is not valid inside --fofn input (%s).", source_name);
            status = status == FC_INPUT_RESOLVE_OK ? FC_INPUT_RESOLVE_FOFN_ERROR : status;
            continue;
        }
        if (fc_stat_path(line, &st) != 0) {
            FC_LOG_ERROR("%s: could not stat '%s': %s", source_name, line, strerror(errno));
            status = status == FC_INPUT_RESOLVE_OK ? FC_INPUT_RESOLVE_IO_ERROR : status;
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            FC_LOG_ERROR("%s: path '%s' is a directory; expected a file path.", source_name, line);
            status = status == FC_INPUT_RESOLVE_OK ? FC_INPUT_RESOLVE_FOFN_ERROR : status;
            continue;
        }
        int rc = append_path(list, line);
        if (rc != 0) {
            FC_LOG_ERROR("Out of memory reading input list.");
            status = FC_INPUT_RESOLVE_IO_ERROR;
            break;
        }
    }

    free(line);
    return status;
}

static fc_input_resolve_status read_fofn_sources(
    const fc_path_list *sources,
    fc_path_list *out_paths)
{
    if (!sources || !sources->paths || sources->n_paths == 0 || !out_paths) {
        return FC_INPUT_RESOLVE_INVALID_ARGUMENT;
    }

    out_paths->paths = NULL;
    out_paths->n_paths = 0;

    fc_input_resolve_status status = FC_INPUT_RESOLVE_OK;
    fc_path_list collected = { 0 };
    bool consumed_stdin = false;

    for (size_t i = 0; i < sources->n_paths; ++i) {
        const char *source = sources->paths[i];
        if (!source || source[0] == '\0') {
            continue;
        }

        if (strcmp(source, "-") == 0) {
            if (consumed_stdin) {
                FC_LOG_ERROR("stdin ('-') can only be listed once with --fofn.");
                status = status == FC_INPUT_RESOLVE_OK ? FC_INPUT_RESOLVE_FOFN_ERROR : status;
                continue;
            }
            consumed_stdin = true;
            fc_input_resolve_status rc = read_paths_from_stream(
                stdin, "stdin", &collected);
            if (rc != FC_INPUT_RESOLVE_OK) {
                status = status == FC_INPUT_RESOLVE_OK ? rc : status;
            }
            continue;
        }

        FILE *fp = fopen(source, "r");
        if (!fp) {
            FC_LOG_ERROR("could not open --fofn source '%s': %s", source, strerror(errno));
            status = status == FC_INPUT_RESOLVE_OK ? FC_INPUT_RESOLVE_FOFN_ERROR : status;
            continue;
        }

        fc_input_resolve_status rc = read_paths_from_stream(
            fp, source, &collected);
        fclose(fp);
        if (rc != FC_INPUT_RESOLVE_OK) {
            status = status == FC_INPUT_RESOLVE_OK ? rc : status;
        }
    }

    if (status != FC_INPUT_RESOLVE_OK) {
        fc_path_list_free(&collected);
        return status;
    }

    *out_paths = collected;
    return 0;
}
