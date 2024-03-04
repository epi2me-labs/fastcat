#include <string.h>
#include <stdbool.h>
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
    meta->rg = "";
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

    // if RG:Z appears in the seq->comment, assume there are SAM tags to parse
    char* res = NULL;
    bool sam_tags = false;
    if (strlen(meta->comment) > 0) {
        res = strstr(meta->comment, "\tRG:Z:");
        if (res != NULL) {
            sam_tags = true;
        }
    }

    char *pch=NULL, *p1=NULL, *p2=NULL;
    char *key=NULL, *keytype=NULL, *value=NULL;

    char sam_token[2] = "\t";
    char fq_token[2] = " ";
    char* token = fq_token;
    if (sam_tags) {
        token = sam_token;
    }
    pch = strtok_r(meta->comment, token, &p1);
    while (pch != NULL) {

        if (sam_tags) {
            // split to tag:type:value
            key = strtok_r(pch, ":", &p2);
            keytype = strtok_r(NULL, ":", &p2);
            value = strtok_r(NULL, "", &p2);
        }
        else {
            // split words on `=`
            key = strtok_r(pch, "=", &p2);
            keytype = NULL;
            value = strtok_r(NULL, "", &p2);
        }

        // if there was no delimiter in the word, value will be NULL --> add word to `rest`
        if (value == NULL && !sam_tags) {
            ksprintf_with_opt_delim(meta->rest, " ", "%s", key);
        } else {
            if (!strcmp(key, "runid") || !strcmp(key, "RD")) {
                // we'll output RD depending on the value of RG, later
                meta->runid = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "RD:Z:%s", meta->runid);
            }
            else if (!strcmp(key, "RG")) {
                meta->rg = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "RG:Z:%s", value);
            }
            else if (!strcmp(key, "flow_cell_id") || !strcmp(key, "FC")) {
                meta->flow_cell_id = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "FC:Z:%s", value);
            }
            else if (!strcmp(key, "barcode") || !strcmp(key, "BC")) {
                meta->barcode = value;
                meta->ibarcode = atoi(value+7);  // "unclassified" -> 0
                ksprintf_with_opt_delim(meta->tags_str, "\t", "BC:Z:%s", value);
            }
            else if (!strcmp(key, "barcode_alias") || !strcmp(key, "BA")) {
                meta->barcode_alias = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "BA:Z:%s", value);
            }
            else if (!strcmp(key, "read") || !strcmp(key, "RN") || !strcmp(key, "rn")) {
                meta->read_number = atoi(value);
                ksprintf_with_opt_delim(meta->tags_str, "\t", "rn:i:%s", value);
            }
            else if (!strcmp(key, "CH") || !strcmp(key, "ch")) {
                meta->channel = atoi(value);
                ksprintf_with_opt_delim(meta->tags_str, "\t", "ch:i:%s", value);
            }
            else if (!strcmp(key, "start_time") || !strcmp(key, "ST") || !strcmp(key, "st")) {
                meta->start_time = value;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "st:Z:%s", value);
            } else {
                if (sam_tags) {
                    // pass through all other tags
                    ksprintf_with_opt_delim(meta->tags_str, "\t", "%s:%s:%s", key, keytype, value);
                }
                else {
                    // long form key=value was not mapped to a SAM tag, send it to CO via meta->rest
                    ksprintf_with_opt_delim(meta->rest, " ", "%s=%s", key, value);
                }
            }
        }
        pch = strtok_r(NULL, token, &p1);
    }

    // if there is a `rest`, add it to the tags (also check that the first char of rest
    // is not ' ', in which case something must have gone wrong)
    if (meta->rest->l != 0 && meta->rest->s[0] != ' ') {
        ksprintf_with_opt_delim(meta->tags_str, "\t", "CO:Z:%s", meta->rest->s);
    }

    // Populate RD (runid) from RG if present, and RD was not tagged.
    if (strlen(meta->runid) == 0) {
        if (strlen(meta->rg) > 0) {
            char* runid = strtok(meta->rg, "_");
            // We'll at least check the runid is the length that we will accept before
            // taking it at face value for inclusion downstream
            if (strlen(runid) == 40) {
                meta->runid = runid;
                ksprintf_with_opt_delim(meta->tags_str, "\t", "RD:Z:%s", runid);
            }
        }
    }

    return meta;
}
