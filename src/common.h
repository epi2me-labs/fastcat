#ifndef _MODBAMBED_COMMON_H
#define _MODBAMBED_COMMON_H

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


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg);


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length);

// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
void kahan_sum(double* sum, double term, double* c);

float mean_qual(char* qual, size_t len);
float mean_qual_from_bam(u_int8_t* qual, size_t len);

#endif
