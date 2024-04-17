#ifndef FASTCAT_ARGS_H
#define FASTCAT_ARGS_H

typedef struct arguments {
    char *perread;
    char *perfile;
    char *runids;
    char *sample;
    size_t min_length;
    size_t max_length;
    float min_qscore;
    int recurse;
    size_t reheader;
    char* demultiplex_dir;
    char* histograms;
    char **files;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
