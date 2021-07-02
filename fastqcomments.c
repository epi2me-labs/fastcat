#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "kstring.h"
#include "fastqcomments.h"

// The caller is responsible for freeing meta and meta->comment
read_meta parse_read_meta(kstring_t comment) {
    char* data = calloc(comment.l, sizeof(char));
    strncpy(data, comment.s, comment.l);
    read_meta meta = calloc(1, sizeof (_read_meta));
    meta->comment = data;
    meta->runid = "";
    meta->flow_cell_id = "";
    meta->barcode = "";
    meta->ibarcode = 0;
    meta->barcode_alias = "";
    meta->read_number = 0;

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
            }
            else if (!strcmp(key, "flow_cell_id")) {
                meta->flow_cell_id = pch;
            }
            else if (!strcmp(key, "barcode")) {
                meta->barcode = pch;
                meta->ibarcode = atoi(pch+7);
            }
            else if (!strcmp(key, "barcode_alias")) {
                meta->barcode_alias = pch;
            }
            else if (!strcmp(key, "read")) {
                meta->read_number = atoi(pch);
            }
            key = NULL;
        }
        pch = strtok(NULL, " =");
    }
    return meta;
}


void destroy_read_meta(read_meta meta) {
    free(meta->comment);
    free(meta);
}
