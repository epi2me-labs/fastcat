#include "core/fc_failures.h"

const char *failure_type[NUM_FAILURE_CODES] = {
#define X(code) #code,
    PARSE_CODES
#undef X
};
