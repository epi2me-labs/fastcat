#ifndef _FASTCAT_REGITER_H
#define _FASTCAT_REGITER_H

#include <stddef.h>

#include "htslib/sam.h"

typedef struct {
    char *chr;
    int tid;
    int64_t start;
    int64_t end;
} _bed_region;

typedef _bed_region *bed_region;

typedef struct {
    _bed_region *regions;
    size_t n_regions;
    size_t _capacity;
} _bed_regions;

typedef _bed_regions *bed_regions;

typedef struct {
    char *chr;
    int start;
    int end;
    int tid;
    bed_regions regions;
    size_t region_index;
    const sam_hdr_t *hdr;
    int error; // &1: couldn't open BED,
               // &2: couldn't parse single_region
} regiter;

typedef enum {
    REGITER_OK = 0,
    REGITER_EOF = -1,
    REGITER_PARSE_ERROR = -2,
    REGITER_REF_MISSING = -3
} regiter_status;

int region_from_string(char *input, char **chr, int *start, int *end);
char *region_to_string(const bed_region region);
int region_from_bed(FILE *bed_fp, char **chr, int *start, int *end);
regiter init_region_iterator(const char *bed_file, const char *single_region, const sam_hdr_t *hdr);
void destroy_region_iterator(regiter *it);
regiter_status next_region(regiter *it);

bed_regions init_bed(const char *bed_path, const sam_hdr_t *hdr);
bed_regions init_bed_from_sam(const sam_hdr_t *hdr, uint32_t segment_length);
void bed_regions_merge_overlaps(bed_regions regions);
void destroy_bed(bed_regions regions);

#endif
