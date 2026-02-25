#ifndef FASTCAT_FC_STATS_ENGINE_H
#define FASTCAT_FC_STATS_ENGINE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fc_read.h"
#include "stats.h"

typedef struct {
    read_stats *length_stats;
    read_stats *qual_stats;
    read_stats *acc_stats;
    read_stats *cov_stats;
    read_stats *length_stats_unmapped;
    read_stats *qual_stats_unmapped;
    read_stats *polya_stats;
    float polya_cover;
    float polya_qual;
    bool polya_rev;
    bool track_alignment;
} fc_stats_engine;

typedef struct {
    uint32_t read_length;
    float mean_quality;
    size_t qstart;
    size_t qend;
    size_t rstart;
    size_t rend;
    size_t aligned_ref_len;
    size_t ref_length;
    size_t match;
    size_t ins;
    size_t delt;
    size_t sub;
    size_t length;
    float coverage;
    float ref_cover;
    float iden;
    float acc;
    char direction;
} fc_stats_aligned;

typedef enum {
    FC_STATS_SKIP = 0,
    FC_STATS_UNMAPPED,
    FC_STATS_ALIGNED
} fc_stats_read_kind;

fc_stats_engine *fc_stats_engine_create(
    bool track_alignment,
    bool track_polya,
    float polya_cover,
    float polya_qual,
    bool polya_rev);
void fc_stats_engine_destroy(fc_stats_engine *engine);

fc_stats_read_kind fc_stats_engine_add_read(
    fc_stats_engine *engine, const fc_read *fr,
    const sam_hdr_t *hdr, fc_stats_aligned *out);

#endif
