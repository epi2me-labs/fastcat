#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "kstring.h"
#include "fastqcomments.h"

// The caller is responsible for calling destroy_read_meta on the returned object.
read_meta parse_read_meta(kstring_t comment) {
    char* data = calloc(comment.l + 1, sizeof(char));
    strncpy(data, comment.s, comment.l);
    read_meta meta = calloc(1, sizeof (_read_meta));
    meta->comment = data;
    meta->runid = "";
    meta->flow_cell_id = "";
    meta->barcode = "";
    meta->ibarcode = 0;
    meta->barcode_alias = "";
    meta->start_time = "";
    meta->read_number = 0;
    meta->channel = 0;
    meta->valid = 0;  // tracks all fields present
    size_t nfields = 7;

    char* pch;
    pch = strtok(meta->comment, " =");
    char* key = NULL;
    while (pch != NULL) {
        if (key == NULL) {
            key = pch;
        }
        else {
            if (!strcmp(key, "runid")) {
                meta->runid = pch;
                meta->valid += 1;
            }
            else if (!strcmp(key, "flow_cell_id")) {
                meta->flow_cell_id = pch;
                meta->valid += 1;
            }
            else if (!strcmp(key, "barcode")) {
                meta->barcode = pch;
                meta->ibarcode = atoi(pch+7);  // "unclassified" -> 0
                meta->valid += 1;
            }
            else if (!strcmp(key, "barcode_alias")) {
                meta->barcode_alias = pch;
                meta->valid += 1;
            }
            else if (!strcmp(key, "read")) {
                meta->read_number = atoi(pch);
                meta->valid += 1;
            }
            else if (!strcmp(key, "ch")) {
                meta->channel = atoi(pch);
                meta->valid += 1;
            }
            else if (!strcmp(key, "start_time")) {
                meta->start_time = pch;
                meta->valid += 1;
            }
            key = NULL;
        }
        pch = strtok(NULL, " =");
    }
    meta->valid = (meta->valid == nfields);
    return meta;
}


void destroy_read_meta(read_meta meta) {
    free(meta->comment);
    free(meta);
}
