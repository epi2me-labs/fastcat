#ifndef _FASTCAT_COMMON_H
#define _FASTCAT_COMMON_H

#include <stdbool.h>
#include <stdint.h>
#include <argp.h>


/** Simple min/max
 * @param a
 * @param b
 *
 * @returns the min/max of a and b
 *
 */
#define min(a, b) ({ \
    typeof (a) _a = (a); \
    typeof (b) _b = (b); \
    _a < _b ? _a : _b; \
})
#define max(a, b) ({ \
    typeof (a) _a = (a); \
    typeof (b) _b = (b); \
    _a > _b ? _a : _b; \
})


// comparison function for qsort of uint32_t
int cmp_u32(const void *a, const void *b);


/** check if a file exists.
 *
 * @param path path to file
 * @returns true if file exists, false otherwise
 *
 */
bool file_exists(const char *path);


/** check if we can make a directory.
 *
 * @param path path to directory
 * @returns  1 if we can make the directory,
 *           0 if it already exists,
 *          -1 if we cannot make it.
 *
 */
int can_make_dir(const char *path);


/** mkdir a directory structure recursively.
 *
 * @param directory path to ensure exists
 *
 */
int mkdir_p(const char* path);


/** mkdir a directory structure recursively, but fail if pre-exists.
 *
 * @param path directory path
 *
 */
int mkdir_hier(char* path);


/** Ensure parent directory exists for a given path.
 *
 *  @param path path to ensure parent directory exists.
 *  @returns 0 on success, -1 on failure.
 *
 *  This function checks if the parent directory of the given path exists,
 *  and creates it if it does not.
 */
int ensure_parent_dir_exists(const char* path);

/** Next two functions used in argument parses to handle space-separated lists of arguments.
 *
 */
void slurp_args(char ***arr, size_t *n, char *first, struct argp_state *state);
void slurp_ints(uint32_t **arr, size_t *n, char *first, struct argp_state *state);

/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg);


/** Reallocates memory with a message on failure.
 *
 *  @param ptr pointer to realloc.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xrealloc(void *ptr, size_t size, char* msg);


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
void *xrecalloc(void* ptr, size_t orig, size_t num, size_t size, char* msg);

/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, size_t position, size_t length);

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
void kahan_sum(double* sum, double term, double* c);

float mean_qual(char* qual, size_t len);
float mean_qual_naive(char* qual, size_t len);
float mean_qual_from_bam(uint8_t* qual, size_t len);
float mean_qual_from_bam_naive(uint8_t* qual, size_t len);

typedef struct readgroup {
    char* readgroup;
    char* runid;
    char* basecaller;
    char* modcaller;
    char* barcode;
} readgroup;

readgroup* create_rg_info(char* rg);
void destroy_rg_info(readgroup* rg);

#endif
