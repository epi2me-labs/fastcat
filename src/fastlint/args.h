#ifndef _FASTLINT_ARGS_H
#define _FASTLINT_ARGS_H

#include <stdbool.h>

#define MAX_INPUT_FILES 1024

typedef struct arguments {
    const char* fastq[MAX_INPUT_FILES + 1];
    int w;
    int t;
    double max_proportion;
    size_t nfiles;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
