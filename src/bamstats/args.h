#ifndef _MODBAMBED_ARGS_H
#define _MODBAMBED_ARGS_H

#include <stdbool.h>


typedef struct arguments {
    const char** bam;
    char* ref;
    char* region;
    char* read_group;
    char tag_name[2];
    int tag_value;
    int threads;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif