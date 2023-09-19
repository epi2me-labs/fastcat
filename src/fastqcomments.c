#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "fastqcomments.h"


// like `ksprintf()`, but will put the optional delimiter `d` before the added string if
// `s` is not empty (and skip the delimiter otherwise)
#define ksprintf_with_opt_delim(s, d, fmt, ...) \
    ksprintf(s, "%s" fmt, s->l == 0 ? "" : d, __VA_ARGS__)


read_meta create_read_meta(const kstring_t* comment) {
    read_meta meta = xalloc(1, sizeof(_read_meta), "meta");
    meta->comment = xalloc(comment->l + 1, sizeof(char), "meta->comment");
    strncpy(meta->comment, comment->s, comment->l);
    meta->runid = "";
    meta->flow_cell_id = "";
    meta->barcode = "";
    meta->ibarcode = 0;
    meta->barcode_alias = "";
    meta->start_time = "";
    meta->read_number = 0;
    meta->channel = 0;
    meta->rest = xalloc(1, sizeof(kstring_t), "meta->rest");
    ks_initialize(meta->rest);
    meta->tags_str = xalloc(1, sizeof(kstring_t), "meta->tags_str");
    ks_initialize(meta->tags_str);

    return meta;
}

void destroy_read_meta(read_meta meta) {
    free(meta->comment);
    free(meta->rest->s);
    free(meta->rest);
    free(meta->tags_str->s);
    free(meta->tags_str);
    free(meta);
}

// The caller is responsible for calling destroy_read_meta on the returned object.
read_meta parse_read_meta(kstring_t comment) {
    read_meta meta = create_read_meta(&comment);

    char *pch=NULL, *key=NULL, *value=NULL, *p1=NULL, *p2=NULL;
    pch = strtok_r(meta->comment, " ", &p1);
    while (pch != NULL) {
        // split words on `=`
        key = strtok_r(pch, "=", &p2);
        value = strtok_r(NULL, "", &p2);

        // if there was no `=` in the word, value will be NULL --> add word to `rest`
        if (value == NULL) {
            ksprintf_with_opt_delim(meta->rest, " ", "%s", key);
        } else {
            if (!strcmp(key, "runid")) {
                meta->runid = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "RD:Z:%s", value);
            }
            else if (!strcmp(key, "flow_cell_id")) {
                meta->flow_cell_id = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "FC:Z:%s", value);
            }
            else if (!strcmp(key, "barcode")) {
                meta->barcode = value;
                meta->ibarcode = atoi(value+7);  // "unclassified" -> 0
                ksprintf_with_opt_delim(meta->tags_str, "\t", "BC:Z:%s", value);
            }
            else if (!strcmp(key, "barcode_alias")) {
                meta->barcode_alias = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "BA:Z:%s", value);
            }
            else if (!strcmp(key, "read")) {
                meta->read_number = atoi(value);
                ksprintf_with_opt_delim(meta->tags_str, "\t", "RN:i:%s", value);
            }
            else if (!strcmp(key, "ch")) {
                meta->channel = atoi(value);
                ksprintf_with_opt_delim(meta->tags_str, "\t", "CH:i:%s", value);
            }
            else if (!strcmp(key, "start_time")) {
                meta->start_time = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "ST:Z:%s", value);
            } else {
                // the word was an unknown key-value pair --> add `key=val` to rest
                ksprintf_with_opt_delim(meta->rest, " ", "%s=%s", key, value);
            }
        }
        pch = strtok_r(NULL, " ", &p1);
    }

    // if there is a `rest`, add it to the tags (also check that the first char of rest
    // is not ' ', in which case something must have gone wrong)
    if (meta->rest->l != 0 && meta->rest->s[0] != ' ') {
        ksprintf_with_opt_delim(meta->tags_str, "\t", "CO:Z:%s", meta->rest->s);
    }

    return meta;
}
