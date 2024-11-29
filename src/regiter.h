#ifndef _FASTCAT_REGITER_H
#define _FASTCAT_REGITER_H

#include "htslib/sam.h"

typedef struct {
    char *chr;
    int start;
    int end;
    FILE *bed_fp;
    char *single_region;
    int n_regions;
    sam_hdr_t *hdr;
    int mode; // 0: BED file, 1: single region
    int error;  // &1: couldn't open BED,
                // &2: couldn't parse single_region
} regiter;


int region_from_string(char* input, char** chr, int* start, int* end);
int region_from_bed(FILE* bed_fp, char** chr, int* start, int* end);
regiter init_region_iterator(const char *bed_file, const char *single_region, sam_hdr_t *hdr);
void destroy_region_iterator(regiter *it);
int next_region(regiter *it);

#endif
