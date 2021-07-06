#ifndef FASTCAT_ARGS_H
#define FASTCAT_ARGS_H

typedef struct arguments {
    char *perread;
    char *perfile;
    char *sample;
    size_t min_length;
    size_t max_length;
    float min_qscore;
    size_t recurse;
    char* demultiplex_dir;
    char **files;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
