#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"


/* The following two functions were adpated from:
 * https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
 */
static int maybe_mkdir(const char* path, mode_t mode) {
    struct stat st;
    errno = 0;

    // Try to make the directory
    if (mkdir(path, mode) == 0)
        return 0;

    // If it fails for any reason but EEXIST, fail
    if (errno != EEXIST)
        return -1;

    // Check if the existing path is a directory
    if (stat(path, &st) != 0)
        return -1;

    // If not, fail with ENOTDIR
    if (!S_ISDIR(st.st_mode)) {
        errno = ENOTDIR;
        return -1;
    }

    errno = 0;
    return 0;
}

/** mkdir a directory structure recursively, but fail if pre-exists.
 *
 * @param path directory path to ensure exists
 *
 */
int mkdir_hier(char *path) {

    char *_path = NULL;
    char *p; 
    int result = -1;
    mode_t mode = 0700;
    errno = 0;

    // if we can just make the directory, fine. If it exists
    // already then exit
    if (mkdir(path, mode) != 0 && errno == EEXIST)
        return -1;
        
    _path = strdup(path);
    if (_path == NULL)
        goto out;

    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (maybe_mkdir(_path, mode) != 0)
                goto out;
            *p = '/';
        }
    }   

    if (maybe_mkdir(_path, mode) != 0)
        goto out;

    result = 0;
out:
    free(_path);
    return result;
}




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
char *substring(char *string, size_t position, size_t length) {
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

int replace_char(char *str, char orig, char rep) {
    char *ix = str;
    int n = 0;
    while((ix = strchr(ix, orig)) != NULL) {
        *ix++ = rep;
        n++;
    }
    return n;
}

const double qprobs[100] = {
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
    2.51188643e-10, 1.99526231e-10, 1.58489319e-10, 1.25892541e-10};


void kahan_sum(double* sum, double term, double* c) {
    double y = term + *c;
    double t = *sum + y;
    *c = (t - *sum) - y;
    *sum = t;
}


float mean_qual(char* qual, size_t len) {
    if (len == 0 ) return nanf("");
    double qsum = 0;
    double c = 0;
    for (size_t i=0; i<len; ++i) {
        int q = (int)(qual[i]) - 33;
        kahan_sum(&qsum, qprobs[q], &c);
    }
    qsum /= len;
    return -10 * log10(qsum);
}

float mean_qual_from_bam(u_int8_t* qual, size_t len) {
    if (len == 0 || qual[0] == 0xff ) return nanf("");
    double qsum = 0;
    double c = 0;
    for (size_t i=0; i<len; ++i) {
        int q = (int)(qual[i]);
        kahan_sum(&qsum, qprobs[q], &c);
    }
    qsum /= len;
    return -10 * log10(qsum);
}
