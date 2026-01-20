#ifndef FASTCAT_FC_FILTER_H
#define FASTCAT_FC_FILTER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "htslib/sam.h"

#define FC_NOTPRIMARY (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP)
#define FC_PRIMARY_QC_EXCLUDE (FC_NOTPRIMARY | BAM_FQCFAIL | BAM_FDUP)

static inline bool fc_is_unmapped(const bam1_t *b) {
    return (b->core.flag & BAM_FUNMAP) != 0;
}

static inline bool fc_is_primary(const bam1_t *b) {
    return (b->core.flag & FC_NOTPRIMARY) == 0;
}

static inline bool fc_is_primary_qc_pass(const bam1_t *b) {
    return (b->core.flag & FC_PRIMARY_QC_EXCLUDE) == 0;
}

typedef struct {
    size_t min_length;
    size_t max_length;
    float min_qscore;
} fc_filter_opts;

typedef enum {
    FC_FILTER_OK = 0,
    FC_FILTER_TOO_LONG,
    FC_FILTER_TOO_SHORT,
    FC_FILTER_LOW_QUALITY
} fc_filter_result;

fc_filter_result fc_filter_length(size_t length, const fc_filter_opts *opts);
fc_filter_result fc_filter_quality(float mean_q, const fc_filter_opts *opts);
bool fc_filter_dust(const uint8_t *seq, size_t len, int t, int w, double max_fraction);
bool fc_filter_int_tag(const bam1_t *b, const char tag_name[2], int tag_value, bool keep_missing);
bool fc_filter_rg(const bam1_t *b, const char *read_group, bool keep_missing);

#endif
