#ifndef _MODBAMBED_ARGS_H
#define _MODBAMBED_ARGS_H

#include <stdbool.h>


typedef struct arguments {
    const char** bam;
    char* flagstats;
    char* runids;
    char* basecallers;
    char* histograms;
    char *sample;
    char* ref;
    char* region;
    char* read_group;
    char tag_name[2];
    int tag_value;
    int threads;
    bool unmapped;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
