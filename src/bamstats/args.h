#ifndef _BAMSTATS_ARGS_H
#define _BAMSTATS_ARGS_H

#include <stdbool.h>


typedef struct arguments {
    const char* bam;
    char* flagstats;
    char* runids;
    char* basecallers;
    char* histograms;
    bool poly_a;
    float poly_a_cover;
    float poly_a_qual;
    bool poly_a_rev;
    char *sample;
    char* ref;
    char* region;
    char* bed;
    char* read_group;
    char tag_name[2];
    int tag_value;
    int threads;
    bool unmapped;
    bool force_recalc_qual;
    // coverage calculations
    bool coverage;
    char* coverages;
    char** coverage_beds;
    size_t n_coverage_beds;
    char** coverage_names;
    size_t n_coverage_names;
    uint32_t* coverage_thresholds;
    size_t n_coverage_thresholds;
    uint32_t* segments;
    size_t n_segments;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);
void destroy_args(arguments_t *args);

#endif
