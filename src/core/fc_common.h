#ifndef _FASTCAT_COMMON_H
#define _FASTCAT_COMMON_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "core/fc_argp.h"
#include "core/fc_log.h"

/** Simple min/max
 * @param a
 * @param b
 *
 * @returns the min/max of a and b
 *
 */
#define min(a, b) ({ \
    typeof(a) _a = (a); \
    typeof(b) _b = (b); \
    _a < _b ? _a : _b; \
})
#define max(a, b) ({ \
    typeof(a) _a = (a); \
    typeof(b) _b = (b); \
    _a > _b ? _a : _b; \
})

#define FC_SNPRINTF_CHECK(buf, size, ...) \
    do { \
        size_t _fc_snprintf_size = (size_t)(size); \
        int _fc_snprintf_n = snprintf((buf), _fc_snprintf_size, __VA_ARGS__); \
        if (_fc_snprintf_n < 0 || (size_t)_fc_snprintf_n >= _fc_snprintf_size) { \
            FC_LOG_ERROR("snprintf failed at %s:%d", __FILE__, __LINE__); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

static inline bool fc_localtime(const time_t *time_value, struct tm *result)
{
#ifdef _WIN32
    return localtime_s(result, time_value) == 0;
#else
    return localtime_r(time_value, result) != NULL;
#endif
}

char *fc_asprintf(const char *msg, const char *fmt, ...);

/** Read one line, growing *line as needed.
 *
 * Mirrors POSIX getline: returns the number of bytes read (including a
 * trailing newline), or -1 on EOF or error. The caller owns *line.
 */
ptrdiff_t fc_getline(char **line, size_t *capacity, FILE *stream);

/** Duplicates a string with a message on failure.
 *
 *  @param s string to duplicate (may be NULL).
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to duplicated string, or NULL if input is NULL.
 *
 */
char *fc_strdup(const char *s, const char *msg);
char *fc_strndup(const char *s, size_t len, const char *msg);
bool fc_find_value_for_key(const char *text, const char *key, const char **value, size_t *value_len);

// Maximum barcode slots: unclassified (0) + barcode0001..barcode1024.
#define FC_MAX_BARCODES 1025

// comparison function for qsort of uint32_t
int cmp_u32(const void *a, const void *b);

// comparison function for qsort of double
int cmp_double(const void *a, const void *b);

/** check if a file exists.
 *
 * @param path path to file
 * @returns true if file exists, false otherwise
 *
 */
bool file_exists(const char *path);

/** mkdir a directory structure recursively.
 *
 * @param directory path to ensure exists
 *
 */
int mkdir_p(const char *path);

/** mkdir a directory structure recursively.
 *
 * @param path directory path
 * @param fail_on_exist if true, return -1 if directory already exists; if false, return 0 if directory exists
 *
 */
int mkdir_hier(char *path, bool fail_on_exist);

/** Remove a file or directory tree. Symbolic links are not followed on POSIX.
 *
 * @param path path to remove; a missing path is treated as success.
 * @returns 0 on success, -1 on failure.
 */
int fc_remove_tree(const char *path);

/** Ensure parent directory exists for a given path.
 *
 *  @param path path to ensure parent directory exists.
 *  @returns 0 on success, -1 on failure.
 *
 *  This function checks if the parent directory of the given path exists,
 *  and creates it if it does not.
 */
int ensure_parent_dir_exists(const char *path);

// Path helpers shared across tools and sinks
char *fc_strip_path(const char *input);
char *fc_join_path(const char *dir, const char *leaf);

/** Next two functions used in argument parses to handle space-separated lists of arguments.
 *
 */
void slurp_args(const char ***arr, size_t *n, const char *first, struct argp_state *state);
void slurp_ints(uint32_t **arr, size_t *n, char *first, struct argp_state *state);

// Strict parsers for user-facing numeric CLI arguments.
bool fc_parse_int_strict(const char *text, int *out);
bool fc_parse_size_strict(const char *text, size_t *out);
bool fc_parse_u32_strict(const char *text, uint32_t *out);
bool fc_parse_float_strict(const char *text, float *out);
bool fc_parse_double_strict(const char *text, double *out);

/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *fc_alloc(size_t num, size_t size, const char *msg);

/** Reallocates memory with a message on failure.
 *
 *  @param ptr pointer to realloc.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *fc_realloc(void *ptr, size_t size, const char *msg);

/** Reallocates memory with a message on failure, zero-initialising the new memory.
 *
 *  @param ptr pointer to realloc.
 *  @param orig original number of elements.
 *  @param num new (total) number of elements.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *fc_recalloc(void *ptr, size_t orig, size_t num, size_t size, const char *msg);

/** Globally replace a char in a char*
 *
 * @param str char* source string
 * @param orig original character
 * @param rep replacement
 * @returns number of times replacement made
 *
 */
int replace_char(char *str, char orig, char rep);

// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
void kahan_sum(double *sum, double term, double *c);

float mean_qual(const char *qual, size_t len);
float mean_qual_naive(const char *qual, size_t len);
float mean_qual_from_bam(const uint8_t *qual, size_t len);
float mean_qual_from_bam_naive(const uint8_t *qual, size_t len);

typedef struct readgroup {
    char *readgroup;
    char *runid;
    char *basecaller;
    char *modcaller;
    char *barcode;
} readgroup;

readgroup *create_rg_info(char *rg);
void destroy_rg_info(readgroup *rg);

#endif
