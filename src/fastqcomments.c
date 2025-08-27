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
    meta->rg_info = NULL;
    meta->runid = "";
    meta->basecaller = "";
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
    destroy_rg_info(meta->rg_info);
    free(meta->rest->s);
    free(meta->rest);
    free(meta->tags_str->s);
    free(meta->tags_str);
    free(meta);
}

// The caller is responsible for calling destroy_read_meta on the returned object.
read_meta parse_read_meta(kstring_t comment) {
    read_meta meta = create_read_meta(&comment);

    // if an RG or RD tag appears in the seq->comment, assume there are SAM tags to parse
    char* res = NULL;
    bool sam_tags = false;
    if (strlen(meta->comment) > 0) {
        // check if comment starts with "RG:Z:" or "RD:Z:"
        if (!strncmp(meta->comment, "RG:Z:", 5) || !strncmp(meta->comment, "RD:Z:", 5)) {
            sam_tags = true;
        }
        // RG or RD could also appear later in the comment (we include '\t' in the check
        // to be extra stringent)
        res = strstr(meta->comment, "\tRG:Z:");
        if (res != NULL) {
            sam_tags = true;
        }
        res = strstr(meta->comment, "\tRD:Z:");
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
            // we allow empty tags (e.g. 'RG:Z:'); in this case, `keytype` will be
            // non-null, but `value` will be null; we set it to ""
            if (keytype != NULL && value == NULL) value = "";
        }
        else {
            // split words on `=`
            key = strtok_r(pch, "=", &p2);
            bool has_delim = p2 != NULL;  // p2 will be NULL if no = was found
            keytype = NULL;
            value = strtok_r(NULL, "", &p2);
            // similarly, allow empty tags (e.g. barcode=) to ensure that k=v
            // records are appropriately formatted in CO:Z even if blank
            if (value == NULL && has_delim) value = "";
        }

        // if there was no delimiter in the word, value will be NULL --> add word to `rest`
        if (value == NULL) {
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
            // CW-4766 - inconsistent naming of basecall model version id by guppy/minknow/dorado
            else if (!strcmp(key, "basecall_model_version_id") || !strcmp(key, "model_version_id")) {
                meta->basecaller = value;
                // there's no discrete tag defined by guppy/minknow/doroado
                // for this; so not added to `tags_str` (but to `rest` instead)
                ksprintf_with_opt_delim(meta->rest, " ", "%s=%s", key, value);
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

    // if there is a `rest`
    // (also check that the first char of rest is not ' ', in which case something must have gone wrong)
    // first replace all tabs with space to avoid all manner of confusion with Martin
    for (size_t i=0; i<meta->rest->l; i++) {
        if (meta->rest->s[i] == '\t') {
            meta->rest->s[i] = ' ';
        }
    }
    if (meta->rest->l != 0 && meta->rest->s[0] != ' ') {
        ksprintf_with_opt_delim(meta->tags_str, "\t", "CO:Z:%s", meta->rest->s);
    }

    bool need_run_id = strlen(meta->runid) == 0;
    bool need_basecaller = strlen(meta->basecaller) == 0;
    if(strlen(meta->rg) > 0 && (need_run_id || need_basecaller)) {
        readgroup* rg_info = create_rg_info(meta->rg);
        if (need_run_id && rg_info->runid != NULL) {
            meta->runid = rg_info->runid;
            ksprintf_with_opt_delim(meta->tags_str, "\t", "RD:Z:%s", rg_info->runid);
        }
        if (need_basecaller && rg_info->basecaller != NULL) {
            meta->basecaller = rg_info->basecaller;
        }
        meta->rg_info = rg_info;
    }

    return meta;
}
