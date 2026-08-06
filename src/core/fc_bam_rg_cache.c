#include "core/fc_bam_rg_cache.h"

#include <string.h>

#include "htslib/kstring.h"

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/kh_counter.h"

struct fc_bam_rg_cache {
    fc_bam_rg_entry *entries;
    size_t n;
    size_t cap;
    kh_counter_t *index; // we're abusing kh_counter_t as a string->int map, gets around clang static analyzer false positives
};

static void fc_bam_rg_entry_destroy(fc_bam_rg_entry *entry)
{
    if (NULL == entry) {
        return;
    }
    free(entry->id);
    free(entry->ds);
    free(entry->runid);
    free(entry->basecall_model);
    free(entry->modbase_models);
    free(entry->dt);
    free(entry->lb);
    free(entry->pu);
    free(entry->sm);
    free(entry->al);
    destroy_rg_info(entry->parsed_id);
    memset(entry, 0, sizeof(*entry));
}

static char *fc_bam_rg_dup_tag(const sam_hdr_t *hdr, int pos, const char *tag, const char *label)
{
    kstring_t ks = { 0, 0, NULL };
    char *value = NULL;
    if (sam_hdr_find_tag_pos((sam_hdr_t *)hdr, "RG", pos, tag, &ks) >= 0 && ks.l > 0) {
        value = fc_strdup(ks.s, label);
    }
    free(ks.s);
    return value;
}

static void fc_bam_rg_cache_grow(fc_bam_rg_cache *cache)
{
    if (cache->n < cache->cap) {
        return;
    }
    size_t next = (0 == cache->cap) ? 8 : cache->cap * 2;
    cache->entries = fc_recalloc(cache->entries, cache->cap, next, sizeof(*cache->entries), "bam rg cache entries");
    cache->cap = next;
}

fc_bam_rg_cache *fc_bam_rg_cache_create(const sam_hdr_t *hdr)
{
    fc_bam_rg_cache *cache = fc_alloc(1, sizeof(*cache), "bam rg cache");
    cache->index = kh_counter_init();
    if (NULL == hdr) {
        return cache;
    }

    int n_rg = sam_hdr_count_lines((sam_hdr_t *)hdr, "RG");
    for (int pos = 0; pos < n_rg; ++pos) {
        const char *id = sam_hdr_line_name((sam_hdr_t *)hdr, "RG", pos);
        if (NULL == id || '\0' == id[0]) {
            continue;
        }

        if (0 != kh_counter_val(cache->index, (char *)id)) {
            FC_LOG_WARN("Duplicate @RG ID '%s' in BAM header, keeping first entry.", id);
            continue;
        }

        fc_bam_rg_cache_grow(cache);
        fc_bam_rg_entry *entry = &cache->entries[cache->n];
        entry->id = fc_strdup(id, "bam rg id");
        entry->ds = fc_bam_rg_dup_tag(hdr, pos, "DS", "bam rg ds");
        entry->dt = fc_bam_rg_dup_tag(hdr, pos, "DT", "bam rg dt");
        entry->lb = fc_bam_rg_dup_tag(hdr, pos, "LB", "bam rg lb");
        entry->pu = fc_bam_rg_dup_tag(hdr, pos, "PU", "bam rg pu");
        entry->sm = fc_bam_rg_dup_tag(hdr, pos, "SM", "bam rg sm");
        entry->al = fc_bam_rg_dup_tag(hdr, pos, "al", "bam rg al");

        if (NULL != entry->ds) {
            const char *value = NULL;
            size_t value_len = 0;
            if (fc_find_value_for_key(entry->ds, "runid", &value, &value_len) && value_len > 0) {
                entry->runid = fc_strndup(value, value_len, "bam rg runid");
            }
            if (fc_find_value_for_key(entry->ds, "basecall_model", &value, &value_len) && value_len > 0) {
                entry->basecall_model = fc_strndup(value, value_len, "bam rg basecall model");
            }
            if (fc_find_value_for_key(entry->ds, "modbase_models", &value, &value_len) && value_len > 0) {
                entry->modbase_models = fc_strndup(value, value_len, "bam rg modbase models");
            }
        }

        if ((size_t)INT_MAX <= cache->n) {
            fc_bam_rg_cache_destroy(cache);
            return NULL;
        }
        if ((size_t)-1 == kh_counter_add(cache->index, entry->id, (int)cache->n + 1)) {
            fc_bam_rg_cache_destroy(cache);
            return NULL;
        }
        cache->n += 1;
    }

    return cache;
}

void fc_bam_rg_cache_destroy(fc_bam_rg_cache *cache)
{
    if (NULL == cache) {
        return;
    }
    for (size_t i = 0; i < cache->n; ++i) {
        fc_bam_rg_entry_destroy(&cache->entries[i]);
    }
    free(cache->entries);
    kh_counter_destroy(cache->index);
    free(cache);
}

const fc_bam_rg_entry *fc_bam_rg_cache_find(const fc_bam_rg_cache *cache, const char *rg_id)
{
    if (NULL == cache || NULL == cache->index || NULL == rg_id || '\0' == rg_id[0]) {
        return NULL;
    }
    int slot = kh_counter_val(cache->index, (char *)rg_id);
    if (0 >= slot) {
        return NULL;
    }
    return &cache->entries[(size_t)slot - 1];
}

const fc_bam_rg_entry *fc_bam_rg_cache_resolve(fc_bam_rg_cache *cache, const char *rg_id)
{
    if (NULL == cache || NULL == cache->index || NULL == rg_id || '\0' == rg_id[0]) {
        return NULL;
    }

    int slot = kh_counter_val(cache->index, (char *)rg_id);
    fc_bam_rg_entry *entry = NULL;
    if (slot > 0) {
        entry = &cache->entries[(size_t)slot - 1];
    } else {
        if ((size_t)INT_MAX <= cache->n) {
            return NULL;
        }
        fc_bam_rg_cache_grow(cache);
        entry = &cache->entries[cache->n];
        entry->id = fc_strdup(rg_id, "bam rg id");
        if ((size_t)-1 == kh_counter_add(cache->index, entry->id, (int)cache->n + 1)) {
            fc_bam_rg_entry_destroy(entry);
            return NULL;
        }
        cache->n += 1;
    }

    if (NULL == entry->parsed_id) {
        entry->parsed_id = create_rg_info(entry->id);
    }
    return entry;
}
