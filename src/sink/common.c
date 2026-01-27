#include "sink/common.h"

#include <stdlib.h>
#include <string.h>

char *fc_sink_strip_path(const char *input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (len > 0 && input[len - 1] == '/') len--;
    char *output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
}
