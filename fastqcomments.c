#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "kstring.h"
#include "fastqcomments.h"

// The caller is responsible for freeing meta.comment
read_meta parse_read_meta(kstring_t comment) {
    read_meta meta = {NULL, "", "", "", 0};
    char* data = calloc(comment.l, sizeof(char));
    strncpy(data, comment.s,comment.l);
    meta.comment = data;
    char* pch;
    pch = strtok(meta.comment, " =");
    char* key = NULL;
    while (pch != NULL) {
        if (key == NULL) {
            key = pch;
        }
        else {
            if (!strcmp(key, "runid")) {
                meta.runid = pch;
            }
            else if (!strcmp(key, "flow_cell_id")) {
                meta.flow_cell_id = pch;
            }
            else if (!strcmp(key, "barcode")) {
                meta.barcode = pch;
                meta.ibarcode = atoi(pch+7);
            }
            key = NULL;
        }
        pch = strtok(NULL, " =");
    }
    return meta;
}
