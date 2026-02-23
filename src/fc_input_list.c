#include "fc_input_list.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "fc_log.h"

static int append_path(fc_path_list *list, const char *path) {
    char *copy = strdup(path);
    if (!copy) return ENOMEM;

    char **tmp = realloc(list->paths, (list->n_paths + 1) * sizeof(*tmp));
    if (!tmp) {
        free(copy);
        return ENOMEM;
    }

    list->paths = tmp;
    list->paths[list->n_paths] = copy;
    list->n_paths++;
    return 0;
}


static int read_paths_from_stream(
    FILE *fp,
    const char *source_name,
    fc_file_kind expected_kind,
    fc_path_list *list) {
    int status = 0;
    char *line = NULL;
    size_t cap = 0;
    ssize_t nchr = 0;

    while ((nchr = getline(&line, &cap, fp)) != -1) {
        struct stat st = {0};
        line[strcspn(line, "\r\n")] = 0;
        if (line[0] == '\0') continue;

        if (strcmp(line, "-") == 0) {
            FC_LOG_ERROR("'-' is not valid inside --fofn input (%s).", source_name);
            status = (status == 0) ? EINVAL : status;
            continue;
        }
        if (stat(line, &st) != 0) {
            FC_LOG_ERROR("%s: could not stat '%s': %s", source_name, line, strerror(errno));
            status = (status == 0) ? errno : status;
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            FC_LOG_ERROR("%s: path '%s' is a directory; expected a file path.", source_name, line);
            status = (status == 0) ? EISDIR : status;
            continue;
        }
        if (!fc_path_matches_kind(line, expected_kind)) {
            FC_LOG_ERROR("Input '%s' does not match required input type for this command.", line);
            status = (status == 0) ? EINVAL : status;
            continue;
        }

        int rc = append_path(list, line);
        if (rc != 0) {
            FC_LOG_ERROR("Out of memory reading input list.");
            status = rc;
            break;
        }
    }

    free(line);
    return status;
}


int fc_read_paths_from_fofn_sources(
    const fc_path_list *sources,
    fc_file_kind expected_kind,
    fc_path_list *out_paths) {
    if (!sources || !sources->paths || sources->n_paths == 0 || !out_paths) return EINVAL;

    out_paths->paths = NULL;
    out_paths->n_paths = 0;

    int status = 0;
    fc_path_list collected = {0};
    bool consumed_stdin = false;

    for (size_t i = 0; i < sources->n_paths; ++i) {
        const char *source = sources->paths[i];
        if (!source || source[0] == '\0') continue;

        if (strcmp(source, "-") == 0) {
            if (consumed_stdin) {
                FC_LOG_ERROR("stdin ('-') can only be listed once with --fofn.");
                status = (status == 0) ? EINVAL : status;
                continue;
            }
            consumed_stdin = true;
            int rc = read_paths_from_stream(stdin, "stdin", expected_kind, &collected);
            if (rc != 0) status = (status == 0) ? rc : status;
            continue;
        }

        FILE *fp = fopen(source, "r");
        if (!fp) {
            FC_LOG_ERROR("could not open --fofn source '%s': %s", source, strerror(errno));
            status = (status == 0) ? errno : status;
            continue;
        }

        int rc = read_paths_from_stream(fp, source, expected_kind, &collected);
        fclose(fp);
        if (rc != 0) status = (status == 0) ? rc : status;
    }

    if (status == 0 && collected.n_paths == 0) status = EINVAL;
    if (status != 0) {
        fc_path_list_free(&collected);
        return status;
    }

    *out_paths = collected;
    return 0;
}
