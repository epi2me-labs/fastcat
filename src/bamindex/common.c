#include <stdlib.h>
#include <string.h>
#include <stdio.h>

char* generate_index_filename(const char* input_bam, const char* input_index) {
    char* out_fn;

    if(input_index != NULL) {
        out_fn = malloc(strlen(input_index));
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_index);
    } else {
        out_fn = calloc(sizeof(char), strlen(input_bam) + 5);
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_bam);
        strcat(out_fn, ".bci\0");
    }
    return out_fn;
}
    