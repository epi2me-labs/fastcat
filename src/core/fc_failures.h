#ifndef FASTCAT_FC_FAILURES_H
#define FASTCAT_FC_FAILURES_H

// F_ - file error
// R_ - record error
// missing quality is tracked as an input-format issue
#define PARSE_CODES     \
    X(F_FILE_OK)        \
    X(F_STREAM_ERROR)   \
    X(F_QUAL_MISSING)   \
    X(F_QUAL_TRUNCATED) \
    X(F_UNKNOWN_ERROR)  \
    X(R_RECORD_SEEN)    \
    X(R_RECORD_OK)      \
    X(R_TOO_LONG)       \
    X(R_TOO_SHORT)      \
    X(R_LOW_QUALITY)    \
    X(R_DUST_MASKED)

typedef enum {
#define X(code) code,
    PARSE_CODES
#undef X
        NUM_FAILURE_CODES
} failure_code;

extern const char *failure_type[NUM_FAILURE_CODES];

#endif
