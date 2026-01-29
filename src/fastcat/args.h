#ifndef FASTCAT_ARGS_H
#define FASTCAT_ARGS_H
#include <stdbool.h>


typedef struct arguments {
    char *sample;
    char *output_dir;
    size_t min_length;
    size_t max_length;
    float min_qscore;
    bool dust;
    double max_dust;
    size_t dust_w;
    size_t dust_t;
    int recurse;
    size_t reheader;
    size_t write_bam;
    bool demultiplex;
    char **files;
    size_t reads_per_file;
    int threads;
    bool verbose;
    bool force_error;
    bool per_read;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
