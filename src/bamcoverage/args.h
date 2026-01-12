#ifndef _BAMCOVERAGE_ARGS_H
#define _BAMCOVERAGE_ARGS_H

#include <stdbool.h>


typedef struct arguments {
    char* bam;
    const char** beds;
    size_t n_beds;
    const char** bed_names;
    size_t n_bed_names;
    int threads;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
