#ifndef FASTCAT_FC_TIME_H
#define FASTCAT_FC_TIME_H

#include <stdbool.h>

// Parse a numeric Unix timestamp or the fixed ISO-8601 timestamp prefix
// YYYY-MM-DDTHH:MM:SS.  ISO suffixes (for example .123Z) are accepted for
// compatibility with the existing timestamp parser and do not affect seconds.
bool fc_time_parse_timestamp(const char *text, double *out);

#endif
