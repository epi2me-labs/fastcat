#ifndef _FASTCAT_COMMON_H
#define _FASTCAT_COMMON_H

#include <stdint.h>


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

/** mkdir a directory structure recursively, but fail if pre-exists.
 *
 * @param path directory path to ensure exists
 *
 */
int mkdir_hier(char* path);

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
