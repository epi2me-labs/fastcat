#ifndef FASTCAT_FC_BAM_RG_CACHE_H
#define FASTCAT_FC_BAM_RG_CACHE_H

#include "htslib/sam.h"

typedef struct {
    char *id;
    char *ds;
    char *runid;
    char *basecall_model;
    char *modbase_models;
    char *dt;
    char *lb;
    char *pu;
    char *sm;
    char *al;
} fc_bam_rg_entry;

typedef struct fc_bam_rg_cache fc_bam_rg_cache;

fc_bam_rg_cache *fc_bam_rg_cache_create(const sam_hdr_t *hdr);
void fc_bam_rg_cache_destroy(fc_bam_rg_cache *cache);
const fc_bam_rg_entry *fc_bam_rg_cache_find(const fc_bam_rg_cache *cache, const char *rg_id);

#endif
