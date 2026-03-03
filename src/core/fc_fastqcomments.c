#include "core/fc_fastqcomments.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/sam.h"

#include "core/fc_common.h"
#include "core/fc_log.h"

// like `ksprintf()`, but will put the optional delimiter `d` before the added string if
// `s` is not empty (and skip the delimiter otherwise)
#define ksprintf_with_opt_delim(s, d, fmt, ...) \
    ksprintf((s), "%s" fmt, (s)->l == 0 ? "" : d, __VA_ARGS__)

static bool comment_has_sam_tags(const char *comment)
{
    if (!comment || *comment == '\0') {
        return false;
    }
    const char *p = comment;
    while (*p) {
        while (*p == ' ' || *p == '\t') {
            p++;
        }
        if (!*p) {
            break;
        }
        if (isalnum((unsigned char)p[0]) && isalnum((unsigned char)p[1])
            && p[2] && p[3] && p[4]
            && p[2] == ':' && p[4] == ':') {
            return true;
        }
        while (*p && *p != ' ' && *p != '\t') {
            p++;
        }
    }
    return false;
}

static bool parse_int64(const char *value, int64_t *out)
{
    if (!value || value[0] == '\0') {
        return false;
    }
    errno = 0;
    char *end = NULL;
    long long v = strtoll(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0') {
        return false;
    }
    *out = v;
    return true;
}

static bool has_gpu_prefix(const char *value)
{
    if (!value) {
        return false;
    }
    if (strlen(value) < 4) {
        return false;
    }
    return (tolower((unsigned char)value[0]) == 'g'
        && tolower((unsigned char)value[1]) == 'p'
        && tolower((unsigned char)value[2]) == 'u'
        && value[3] == ':');
}

static int parse_key_value_comment_to_bam(const kstring_t *comment, bam1_t *b)
{
    if (!comment || !comment->s || comment->l == 0 || !b) {
        return 0;
    }

    char *buf = fc_alloc(comment->l + 1, sizeof(char), "comment buf");
    memcpy(buf, comment->s, comment->l);
    buf[comment->l] = '\0';

    kstring_t rest = { 0, 0, NULL };
    ks_initialize(&rest);
    char *runid = NULL;
    char *basecall_model = NULL;
    char *barcode = NULL;
    bool rg_seen = bam_aux_get(b, "RG") != NULL;

    char *pch = NULL, *p1 = NULL;
    pch = strtok_r(buf, " ", &p1);
    while (pch != NULL) {
        // Split on the first '=' without relying on strtok_r("") (UB).
        char *key = pch;
        char *value = NULL;
        char *delim = strchr(pch, '=');
        if (delim != NULL) {
            *delim = '\0';
            value = delim + 1; // May be an empty string.
        }

        if (value == NULL) {
            ksprintf_with_opt_delim(&rest, " ", "%s", key);
        } else if (!strcmp(key, "runid")) {
            runid = value;
        } else if (!strcmp(key, "basecall_model_version_id") || !strcmp(key, "model_version_id")) {
            basecall_model = value;
        } else if (!strcmp(key, "flow_cell_id")) {
            if (bam_aux_update_str(b, "PU", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "sample_id")) {
            if (bam_aux_update_str(b, "LB", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "barcode")) {
            barcode = value;
            if (bam_aux_update_str(b, "SM", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "barcode_alias")) {
            if (bam_aux_update_str(b, "al", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "read")) {
            int64_t ival = 0;
            if (!parse_int64(value, &ival) || ival < 0 || ival > UINT32_MAX) {
                ksprintf_with_opt_delim(&rest, " ", "%s=%s", key, value);
            } else if (bam_aux_update_int(b, "rn", ival) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "ch")) {
            int64_t ival = 0;
            if (!parse_int64(value, &ival) || ival < 0 || ival > UINT32_MAX) {
                ksprintf_with_opt_delim(&rest, " ", "%s=%s", key, value);
            } else if (bam_aux_update_int(b, "ch", ival) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "start_time")) {
            if (bam_aux_update_str(b, "st", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "parent_read_id")) {
            if (bam_aux_update_str(b, "pi", (int)strlen(value), value) < 0) {
                goto err;
            }
        } else if (!strcmp(key, "basecall_gpu")) {
            kstring_t ds = { 0, 0, NULL };
            ks_initialize(&ds);
            if (value[0] != '\0' && !has_gpu_prefix(value)) {
                ksprintf(&ds, "gpu:%s", value);
            } else {
                ksprintf(&ds, "%s", value);
            }
            if (bam_aux_update_str(b, "DS", (int)ds.l, ds.s) < 0) {
                free(ds.s);
                goto err;
            }
            free(ds.s);
        } else {
            ksprintf_with_opt_delim(&rest, " ", "%s=%s", key, value);
        }

        pch = strtok_r(NULL, " ", &p1);
    }

    for (size_t i = 0; i < rest.l; i++) {
        if (rest.s[i] == '\t') {
            rest.s[i] = ' ';
        }
    }
    if (rest.l != 0 && rest.s[0] != ' ') {
        if (bam_aux_update_str(b, "CO", (int)rest.l, rest.s) < 0) {
            goto err;
        }
    }

    // Build an RG tag from runid/model/barcode when no RG is present.
    if (!rg_seen && runid && runid[0] != '\0') {
        kstring_t rg = { 0, 0, NULL };
        ks_initialize(&rg);
        if (basecall_model && basecall_model[0] != '\0') {
            ksprintf(&rg, "%s_%s", runid, basecall_model);
        } else {
            ksprintf(&rg, "%s_unknown", runid);
        }
        if (barcode && barcode[0] != '\0' && strcmp(barcode, "unclassified") != 0) {
            kputc_('_', &rg);
            kputs(barcode, &rg);
        }
        if (rg.l > 0 && bam_aux_update_str(b, "RG", (int)rg.l, rg.s) < 0) {
            free(rg.s);
            goto err;
        }
        free(rg.s);
    }

    free(rest.s);
    free(buf);
    return 0;

err:
    free(rest.s);
    free(buf);
    return -1;
}

int fastq_comment_to_bam(const kstring_t *comment, bam1_t *b)
{
    if (!comment || !comment->s || comment->l == 0 || !b) {
        return 0;
    }

    bool sam_like = comment_has_sam_tags(comment->s);
    if (sam_like) {
        int orig_len = b->l_data;
        if (sam_aux_parse(comment->s, comment->s + comment->l, b, 1) < 0) {
            FC_LOG_WARN("Error parsing SAM-style tags from FASTQ comment.");
            return -1;
        }
        if (b->l_data > orig_len) {
            return 0;
        }
    }

    if (parse_key_value_comment_to_bam(comment, b) < 0) {
        FC_LOG_WARN("Error parsing key=value tags from FASTQ comment.");
        return -1;
    }
    return 0;
}
