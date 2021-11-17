#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg){
    void *res = calloc(num, size);
    if (res == NULL){
        fprintf(stderr, "Failed to allocate mem for %s\n", msg);
        exit(1);
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
void *xrealloc(void *ptr, size_t size, char* msg){
    void *res = realloc(ptr, size);
    if (res == NULL){
        fprintf(stderr, "Failed to reallocate mem for %s\n", msg);
        exit(1);
    }
    return res;
}


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length) {
   char *ptr;
   size_t i;

   ptr = malloc(length + 1);

   for (i = 0 ; i < length ; i++) {
      *(ptr + i) = *(string + position);
      string++;
   }

   *(ptr + i) = '\0';
   return ptr;
}


