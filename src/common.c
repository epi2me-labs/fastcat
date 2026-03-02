#define _GNU_SOURCE
#include "common.h"

#include <assert.h>
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <regex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fc_log.h"

// used with qsort to sort an array of uint32_t
int cmp_u32(const void *a, const void *b)
{
    uint32_t va = *(const uint32_t *)a;
    uint32_t vb = *(const uint32_t *)b;
    if (va < vb) {
        return -1;
    }
    if (va > vb) {
        return 1;
    }
    return 0;
}

/** check if a file exists.
 *
 * @param path path to file
 * @returns true if file exists, false otherwise
 *
 */
bool file_exists(const char *path)
{
    assert(path != NULL);
    struct stat st;
    return stat(path, &st) == 0;
}

/* The following three functions were adpated from:
 * https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
 */
static int maybe_mkdir(const char *path, mode_t mode)
{
    struct stat st;
    errno = 0;

    // Try to make the directory
    if (mkdir(path, mode) == 0) {
        return 0;
    }

    // If it fails for any reason but EEXIST, fail
    if (errno != EEXIST) {
        return -1;
    }

    // Check if the existing path is a directory
    if (stat(path, &st) != 0) {
        return -1;
    }

    // If not, fail with ENOTDIR
    if (!S_ISDIR(st.st_mode)) {
        errno = ENOTDIR;
        return -1;
    }

    errno = 0;
    return 0;
}

/** mkdir a directory structure recursively.
 *
 * @param path directory path to ensure exists
 *
 */
int mkdir_p(const char *path)
{
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */
    char *_path = NULL;
    char *p;
    int result = -1;
    mode_t mode = 0755;
    errno = 0;

    _path = fc_strdup(path, "path copy for mkdir_hier");

    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (maybe_mkdir(_path, mode) != 0) {
                goto out;
            }
            *p = '/';
        }
    }

    if (maybe_mkdir(_path, mode) != 0) {
        goto out;
    }

    result = 0;

out:
    free(_path);
    return result;
}

/** mkdir a directory structure recursively, but fail if pre-exists.
 *
 * @param path directory path
 *
 */
int mkdir_hier(char *path, bool fail_on_exist)
{

    char *_path = NULL;
    char *p;
    int result = -1;
    mode_t mode = 0755;
    errno = 0;

    // if we can just make the directory, fine. If it exists
    // already then fail or succeed depending on fail_on_exist flag
    int first_mkdir_result = mkdir(path, mode);
    int saved_errno = errno;
    if (first_mkdir_result != 0 && saved_errno == EEXIST) {
        return fail_on_exist ? -1 : 0;
    }
    // If first mkdir succeeded, we're done
    if (first_mkdir_result == 0) {
        return 0;
    }

    _path = fc_strdup(path, "path copy for mkdir_p");

    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (maybe_mkdir(_path, mode) != 0) {
                goto out;
            }
            *p = '/';
        }
    }

    if (maybe_mkdir(_path, mode) != 0) {
        goto out;
    }

    result = 0;
out:
    free(_path);
    return result;
}

/** Ensure parent directory exists for a given path.
 *
 *  @param path path to ensure parent directory exists.
 *  @returns 0 on success, -1 on failure.
 *
 *  This function checks if the parent directory of the given path exists,
 *  and creates it if it does not.
 */
int ensure_parent_dir_exists(const char *filepath)
{
    char *path_copy = fc_strdup(filepath, "filepath copy for parent dir check");
    char *dir = dirname(path_copy);
    int rtn = mkdir_p(dir);
    free(path_copy);
    return rtn;
}

static void push_str(const char ***arr, size_t *n, const char *s)
{
    // Use bare realloc + argp_failure to integrate with argp error reporting,
    // rather than fc_realloc which would use FC_LOG_FATAL.
    const char **tmp = realloc(*arr, (*n + 1) * sizeof(**arr));
    if (!tmp) {
        argp_failure(NULL, 1, 0, "out of memory");
        exit(EXIT_FAILURE); // help static analyzer understand control flow
    }
    *arr = tmp;
    (*arr)[(*n)++] = s;
}

bool fc_parse_int_strict(const char *text, int *out)
{
    if (!text || !*text || !out) {
        return false;
    }
    errno = 0;
    char *endp = NULL;
    long val = strtol(text, &endp, 10);
    if (errno == ERANGE || endp == text || *endp != '\0') {
        return false;
    }
    if (val < INT_MIN || val > INT_MAX) {
        return false;
    }
    *out = (int)val;
    return true;
}

bool fc_parse_size_strict(const char *text, size_t *out)
{
    if (!text || !*text || !out) {
        return false;
    }
    if (text[0] == '-') {
        return false;
    }

    errno = 0;
    char *endp = NULL;
    unsigned long long val = strtoull(text, &endp, 10);
    if (errno == ERANGE || endp == text || *endp != '\0') {
        return false;
    }
    if (val > SIZE_MAX) {
        return false;
    }
    *out = (size_t)val;
    return true;
}

bool fc_parse_u32_strict(const char *text, uint32_t *out)
{
    if (!text || !*text || !out) {
        return false;
    }
    if (text[0] == '-') {
        return false;
    }

    errno = 0;
    char *endp = NULL;
    unsigned long long val = strtoull(text, &endp, 10);
    if (errno == ERANGE || endp == text || *endp != '\0') {
        return false;
    }
    if (val > UINT32_MAX) {
        return false;
    }
    *out = (uint32_t)val;
    return true;
}

bool fc_parse_float_strict(const char *text, float *out)
{
    if (!text || !*text || !out) {
        return false;
    }
    errno = 0;
    char *endp = NULL;
    float val = strtof(text, &endp);
    if (errno == ERANGE || endp == text || *endp != '\0') {
        return false;
    }
    if (!isfinite(val)) {
        return false;
    }
    *out = val;
    return true;
}

bool fc_parse_double_strict(const char *text, double *out)
{
    if (!text || !*text || !out) {
        return false;
    }
    errno = 0;
    char *endp = NULL;
    double val = strtod(text, &endp);
    if (errno == ERANGE || endp == text || *endp != '\0') {
        return false;
    }
    if (!isfinite(val)) {
        return false;
    }
    *out = val;
    return true;
}

char *fc_strip_path(const char *input)
{
    if (input == NULL) {
        return NULL;
    }
    size_t len = strlen(input);
    if (len > 0 && input[len - 1] == '/') {
        len--;
    }
    char *output = fc_alloc(len + 1, sizeof(char), "stripped path");
    memcpy(output, input, len);
    output[len] = '\0'; // explicit null termination
    return output;
}

char *fc_join_path(const char *dir, const char *leaf)
{
    if (!dir || !leaf) {
        FC_LOG_FATAL("Cannot join path with NULL components.");
    }
    size_t len = strlen(dir);
    size_t llen = strlen(leaf);
    bool need_sep = len > 0 && dir[len - 1] != '/';

    // Check for overflow before computing total
    if (len > PATH_MAX || llen > PATH_MAX || len + need_sep > PATH_MAX - llen - 1) {
        FC_LOG_FATAL("Path too long: '%s' + '%s' exceeds PATH_MAX (%d)", dir, leaf, PATH_MAX);
    }

    size_t total = len + need_sep + llen + 1;
    char *out = fc_alloc(total, 1, "joined path");
    FC_SNPRINTF_CHECK(out, total, "%s%s%s", dir, need_sep ? "/" : "", leaf);
    return out;
}

void slurp_args(const char ***arr, size_t *n, const char *first, struct argp_state *state)
{
    // slurp up all arguments until the next option
    push_str(arr, n, first);
    while (state->next < state->argc) {
        const char *next = state->argv[state->next];
        if (next[0] == '-') {
            break;
        }
        push_str(arr, n, next);
        state->next++;
    }
}

void slurp_ints(uint32_t **arr, size_t *n, char *first, struct argp_state *state)
{
    uint32_t val = 0;
    if (!fc_parse_u32_strict(first, &val)) {
        argp_error(state, "Invalid integer threshold: '%s'", first);
    }
    // Use bare realloc + argp_failure to integrate with argp error reporting,
    // rather than fc_realloc which would use FC_LOG_FATAL.
    uint32_t *tmp = realloc(*arr, (*n + 1) * sizeof(**arr));
    if (!tmp) {
        argp_failure(state, 1, 0, "out of memory");
        exit(EXIT_FAILURE); // help static analyzer understand control flow
    }
    *arr = tmp;
    (*arr)[(*n)++] = val;

    while (state->next < state->argc) {
        char *next = state->argv[state->next];
        if (next[0] == '-') {
            break;
        }
        if (!fc_parse_u32_strict(next, &val)) {
            argp_error(state, "Invalid integer threshold: '%s'", next);
        }
        // Use bare realloc + argp_failure to integrate with argp error reporting,
        // rather than fc_realloc which would use FC_LOG_FATAL.
        tmp = realloc(*arr, (*n + 1) * sizeof(**arr));
        if (!tmp) {
            argp_failure(state, 1, 0, "out of memory");
            exit(EXIT_FAILURE); // help static analyzer understand control flow
        }
        *arr = tmp;
        (*arr)[(*n)++] = val;
        state->next++;
    }
}

/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *fc_alloc(size_t num, size_t size, const char *msg)
{
    if (num == 0 || size == 0) {
        FC_LOG_FATAL("Invalid zero-size allocation for %s (num=%zu, size=%zu)", msg, num, size);
    }
    if (size > SIZE_MAX / num) {
        FC_LOG_FATAL("Overflow in allocation size for %s (num=%zu, size=%zu)", msg, num, size);
    }
    void *res = calloc(num, size);
    if (res == NULL) {
        FC_LOG_FATAL("Failed to allocate mem for %s", msg);
    }
    return res;
}

/** Reallocates memory with a message on failure.
 *
 *  @param ptr pointer to realloc.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *fc_realloc(void *ptr, size_t size, const char *msg)
{
    void *res = realloc(ptr, size);
    if (res == NULL) {
        FC_LOG_FATAL("Failed to reallocate mem for %s", msg);
    }
    return res;
}

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
void *fc_recalloc(void *ptr, size_t orig, size_t num, size_t size, const char *msg)
{
    if (num < orig) {
        FC_LOG_FATAL("fc_recalloc: new size (%zu) cannot be less than original (%zu) for %s", num, orig, msg);
    }
    if (num == 0 || size == 0) {
        FC_LOG_FATAL("Invalid zero-size reallocation for %s (num=%zu, size=%zu)", msg, num, size);
    }
    if (size > SIZE_MAX / num) {
        FC_LOG_FATAL("Overflow in reallocation size for %s (num=%zu, size=%zu)", msg, num, size);
    }
    if (orig != 0 && size > SIZE_MAX / orig) {
        FC_LOG_FATAL("Overflow in original size calculation for %s (orig=%zu, size=%zu)", msg, orig, size);
    }
    void *res = realloc(ptr, num * size);
    if (res == NULL) {
        FC_LOG_FATAL("Failed to reallocate mem for %s", msg);
    }
    // zero out the new memory (overflow already checked above)
    memset((char *)res + orig * size, 0, (num - orig) * size);
    return res;
}

char *fc_asprintf(const char *msg, const char *fmt, ...)
{
    char *out = NULL;
    va_list ap;
    va_start(ap, fmt);
    int wrote = vasprintf(&out, fmt, ap);
    va_end(ap);
    if (wrote < 0 || out == NULL) {
        FC_LOG_FATAL("asprintf failed for %s", msg ? msg : "formatted string");
    }
    return out;
}

/** Duplicates a string with a message on failure.
 *
 *  @param s string to duplicate (may be NULL).
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to duplicated string, or NULL if input is NULL.
 *
 */
char *fc_strdup(const char *s, const char *msg)
{
    if (!s) {
        return NULL;
    }
    size_t len = strlen(s);
    char *copy = fc_alloc(len + 1, sizeof(char), msg);
    memcpy(copy, s, len + 1);
    return copy;
}

int replace_char(char *str, char orig, char rep)
{
    char *ix = str;
    int n = 0;
    while ((ix = strchr(ix, orig)) != NULL) {
        *ix++ = rep;
        n++;
    }
    return n;
}

// Error probability lookup table: qprobs[q] = 10^(-q/10)
// Extended to 256 entries to cover all possible uint8_t quality scores
const double qprobs[256] = {
    1.00000000e+00, 7.94328235e-01, 6.30957344e-01, 5.01187234e-01,
    3.98107171e-01, 3.16227766e-01, 2.51188643e-01, 1.99526231e-01,
    1.58489319e-01, 1.25892541e-01, 1.00000000e-01, 7.94328235e-02,
    6.30957344e-02, 5.01187234e-02, 3.98107171e-02, 3.16227766e-02,
    2.51188643e-02, 1.99526231e-02, 1.58489319e-02, 1.25892541e-02,
    1.00000000e-02, 7.94328235e-03, 6.30957344e-03, 5.01187234e-03,
    3.98107171e-03, 3.16227766e-03, 2.51188643e-03, 1.99526231e-03,
    1.58489319e-03, 1.25892541e-03, 1.00000000e-03, 7.94328235e-04,
    6.30957344e-04, 5.01187234e-04, 3.98107171e-04, 3.16227766e-04,
    2.51188643e-04, 1.99526231e-04, 1.58489319e-04, 1.25892541e-04,
    1.00000000e-04, 7.94328235e-05, 6.30957344e-05, 5.01187234e-05,
    3.98107171e-05, 3.16227766e-05, 2.51188643e-05, 1.99526231e-05,
    1.58489319e-05, 1.25892541e-05, 1.00000000e-05, 7.94328235e-06,
    6.30957344e-06, 5.01187234e-06, 3.98107171e-06, 3.16227766e-06,
    2.51188643e-06, 1.99526231e-06, 1.58489319e-06, 1.25892541e-06,
    1.00000000e-06, 7.94328235e-07, 6.30957344e-07, 5.01187234e-07,
    3.98107171e-07, 3.16227766e-07, 2.51188643e-07, 1.99526231e-07,
    1.58489319e-07, 1.25892541e-07, 1.00000000e-07, 7.94328235e-08,
    6.30957344e-08, 5.01187234e-08, 3.98107171e-08, 3.16227766e-08,
    2.51188643e-08, 1.99526231e-08, 1.58489319e-08, 1.25892541e-08,
    1.00000000e-08, 7.94328235e-09, 6.30957344e-09, 5.01187234e-09,
    3.98107171e-09, 3.16227766e-09, 2.51188643e-09, 1.99526231e-09,
    1.58489319e-09, 1.25892541e-09, 1.00000000e-09, 7.94328235e-10,
    6.30957344e-10, 5.01187234e-10, 3.98107171e-10, 3.16227766e-10,
    2.51188643e-10, 1.99526231e-10, 1.58489319e-10, 1.25892541e-10,
    1.00000000e-10, 7.94328235e-11, 6.30957344e-11, 5.01187234e-11,
    3.98107171e-11, 3.16227766e-11, 2.51188643e-11, 1.99526231e-11,
    1.58489319e-11, 1.25892541e-11, 1.00000000e-11, 7.94328235e-12,
    6.30957344e-12, 5.01187234e-12, 3.98107171e-12, 3.16227766e-12,
    2.51188643e-12, 1.99526231e-12, 1.58489319e-12, 1.25892541e-12,
    1.00000000e-12, 7.94328235e-13, 6.30957344e-13, 5.01187234e-13,
    3.98107171e-13, 3.16227766e-13, 2.51188643e-13, 1.99526231e-13,
    1.58489319e-13, 1.25892541e-13, 1.00000000e-13, 7.94328235e-14,
    6.30957344e-14, 5.01187234e-14, 3.98107171e-14, 3.16227766e-14,
    2.51188643e-14, 1.99526231e-14, 1.58489319e-14, 1.25892541e-14,
    1.00000000e-14, 7.94328235e-15, 6.30957344e-15, 5.01187234e-15,
    3.98107171e-15, 3.16227766e-15, 2.51188643e-15, 1.99526231e-15,
    1.58489319e-15, 1.25892541e-15, 1.00000000e-15, 7.94328235e-16,
    6.30957344e-16, 5.01187234e-16, 3.98107171e-16, 3.16227766e-16,
    2.51188643e-16, 1.99526231e-16, 1.58489319e-16, 1.25892541e-16,
    1.00000000e-16, 7.94328235e-17, 6.30957344e-17, 5.01187234e-17,
    3.98107171e-17, 3.16227766e-17, 2.51188643e-17, 1.99526231e-17,
    1.58489319e-17, 1.25892541e-17, 1.00000000e-17, 7.94328235e-18,
    6.30957344e-18, 5.01187234e-18, 3.98107171e-18, 3.16227766e-18,
    2.51188643e-18, 1.99526231e-18, 1.58489319e-18, 1.25892541e-18,
    1.00000000e-18, 7.94328235e-19, 6.30957344e-19, 5.01187234e-19,
    3.98107171e-19, 3.16227766e-19, 2.51188643e-19, 1.99526231e-19,
    1.58489319e-19, 1.25892541e-19, 1.00000000e-19, 7.94328235e-20,
    6.30957344e-20, 5.01187234e-20, 3.98107171e-20, 3.16227766e-20,
    2.51188643e-20, 1.99526231e-20, 1.58489319e-20, 1.25892541e-20,
    1.00000000e-20, 7.94328235e-21, 6.30957344e-21, 5.01187234e-21,
    3.98107171e-21, 3.16227766e-21, 2.51188643e-21, 1.99526231e-21,
    1.58489319e-21, 1.25892541e-21, 1.00000000e-21, 7.94328235e-22,
    6.30957344e-22, 5.01187234e-22, 3.98107171e-22, 3.16227766e-22,
    2.51188643e-22, 1.99526231e-22, 1.58489319e-22, 1.25892541e-22,
    1.00000000e-22, 7.94328235e-23, 6.30957344e-23, 5.01187234e-23,
    3.98107171e-23, 3.16227766e-23, 2.51188643e-23, 1.99526231e-23,
    1.58489319e-23, 1.25892541e-23, 1.00000000e-23, 7.94328235e-24,
    6.30957344e-24, 5.01187234e-24, 3.98107171e-24, 3.16227766e-24,
    2.51188643e-24, 1.99526231e-24, 1.58489319e-24, 1.25892541e-24,
    1.00000000e-24, 7.94328235e-25, 6.30957344e-25, 5.01187234e-25,
    3.98107171e-25, 3.16227766e-25, 2.51188643e-25, 1.99526231e-25,
    1.58489319e-25, 1.25892541e-25, 1.00000000e-25, 7.94328235e-26,
    6.30957344e-26, 5.01187234e-26, 3.98107171e-26, 3.16227766e-26
};

inline void kahan_sum(double *sum, double term, double *c)
{
    double y = term - *c;
    double t = *sum + y;
    *c = (t - *sum) - y;
    *sum = t;
}

inline float mean_qual(const char *qual, size_t len)
{
    if (len == 0) {
        return nanf("");
    }
    double qsum = 0;
    double c = 0;
    for (size_t i = 0; i < len; ++i) {
        if ((unsigned char)qual[i] < 33) {
            FC_LOG_FATAL("Invalid FASTQ quality character %d (< 33) at position %zu",
                (unsigned char)qual[i], i);
        }
        uint8_t q = (uint8_t)qual[i] - 33;
        kahan_sum(&qsum, qprobs[q], &c);
    }
    qsum /= len;
    return -10 * log10(qsum);
}

inline float mean_qual_from_bam(const uint8_t *qual, size_t len)
{
    if (len == 0 || qual[0] == 0xff) {
        return nanf("");
    }
    double qsum = 0;
    double c = 0;
    for (size_t i = 0; i < len; ++i) {
        kahan_sum(&qsum, qprobs[qual[i]], &c);
    }
    qsum /= len;
    return -10 * log10(qsum);
}

inline float mean_qual_naive(const char *qual, size_t len)
{
    if (len == 0) {
        return nanf("");
    }
    double qsum = 0;
    for (size_t i = 0; i < len; ++i) {
        if ((unsigned char)qual[i] < 33) {
            FC_LOG_FATAL("Invalid FASTQ quality character %d (< 33) at position %zu",
                (unsigned char)qual[i], i);
        }
        uint8_t q = (uint8_t)qual[i] - 33;
        qsum += qprobs[q];
    }
    qsum /= len;
    return -10 * log10(qsum);
}

inline float mean_qual_from_bam_naive(const uint8_t *qual, size_t len)
{
    if (len == 0 || qual[0] == 0xff) {
        return nanf("");
    }
    double qsum = 0;
    for (size_t i = 0; i < len; ++i) {
        qsum += qprobs[qual[i]];
    }
    qsum /= len;
    return -10 * log10(qsum);
}

// Strip hexadecimal suffixes that samtools merge can add to RG IDs
static void strip_hex_suffix(char *str)
{
    regex_t regex;
    regmatch_t matches[1];

    // samtools formats this as "%s-%0lX" which is a long formatted
    // as hex with leading zeros. BUT there's no field width specified!
    // https://github.com/samtools/samtools/issues/2086
    // This will strip even a solitary `-` from the end of the string,
    // we could be probabilistic and require say 2 hex digits leaving
    // the 16 edge cases of 0 ("") to F. No one should have and RG ID
    // ending in "-", right?
    if (regcomp(&regex, "-[0-9A-Fa-f]{0,8}$", REG_EXTENDED) != 0) {
        FC_LOG_FATAL("Failed to compile regex for RG ID suffix stripping (internal error).");
    }

    // Check if the string matches the pattern, and strip the match from the src
    if (regexec(&regex, str, 1, matches, 0) == 0) {
        str[matches[0].rm_so] = '\0';
    }

    regfree(&regex);
}

void destroy_rg_info(readgroup *rg)
{
    if (rg != NULL) {
        free(rg->readgroup);
        free(rg);
    }
}

// rg is of the form:
// <runid>_<basecalling_model>_<barcode_arrangement>
//
// where:
//   - runid is either (see CW-4704):
//        - a 40 character string representing an acquisition_id sha
//        - a 36 character string representing a protocol_run_id uuid
//   - basecalling_model is a string maybe containing `_`, and containing one or more `@`
//   -    mod_caller is optional part of this starting with `_` after the first `@`
//   - barcode_arrangement is a optional(!) string with an unknown format, but hopefully no `@`
//
// The function always returns an object with a copy of the input. The subfields may be
// NULL pointers if the parsing was incomplete. Client code should therefore always check that
// members are not NULL before use.
//
readgroup *create_rg_info(char *rg)
{
    readgroup *rg_info = fc_alloc(1, sizeof(readgroup), "readgroup");
    rg_info->readgroup = fc_strdup(rg, "readgroup string");
    rg_info->runid = NULL;
    rg_info->basecaller = NULL;
    rg_info->modcaller = NULL;
    rg_info->barcode = NULL;

    // first strip of `-ABCDEF` from the end
    strip_hex_suffix(rg_info->readgroup);

    // I tried to do this with regex, but even chatGPT couldn't give me a
    // POSIX regex that would work. So we'll do it manually (and somewhat
    // more understandably/controllably)

    // runid runs to first `_`
    rg_info->runid = rg_info->readgroup;
    char *delim = strchr(rg_info->runid, '_');
    if (delim == NULL) {
        return rg_info;
    }
    delim[0] = '\0';
    // ensure runid is long enough to be an acquisition sha or protocol uuid
    int runid_l = strlen(rg_info->runid);
    if (runid_l != 36 && runid_l != 40) {
        // free the mutated copy, and reset
        free(rg_info->readgroup);
        rg_info->readgroup = fc_strdup(rg, "readgroup string (reset)");
        rg_info->runid = NULL;
        return rg_info;
    }
    // basecaller + modcaller runs to first `_` after last `@`
    // though barcode is optional, so there may not
    // be a `_` after the last `@`
    rg_info->basecaller = delim + 1;
    delim = strrchr(rg_info->basecaller, '@');
    if (delim == NULL) {
        rg_info->basecaller = NULL;
        return rg_info;
    }
    // modcaller is optional, we can detect its presence by more than one `@` in basecaller
    char *delim1 = strchr(rg_info->basecaller, '@');
    if (delim1 == delim) { // only one `@`
        rg_info->modcaller = NULL;
    } else {
        // if RG string contains an extra '@' separators without a modcaller "_"
        // delimiter after the first '@'. Treat those as "no modcaller" instead of
        // dereferencing a NULL separator.
        char *mod_sep = strchr(delim1, '_'); // modcaller starts at `_` after first `@`
        if (mod_sep != NULL && mod_sep < delim) {
            mod_sep[0] = '\0';
            rg_info->modcaller = mod_sep + 1;
        } else {
            rg_info->modcaller = NULL;
        }
    }
    delim = strchr(delim, '_');
    // barcode is optional
    if (delim) {
        delim[0] = '\0';
        rg_info->barcode = delim + 1;
    }
    return rg_info;
}
