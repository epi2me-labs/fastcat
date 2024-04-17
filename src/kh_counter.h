#ifndef _KHCOUNTER_H
#define _KHCOUNTER_H

#include "htslib/khash.h"


KHASH_MAP_INIT_STR(KH_COUNTER, int)
#define kh_counter_t khash_t(KH_COUNTER)

// create a counter
kh_counter_t *kh_counter_init(void);

// Get a value from a counter 
int kh_counter_val(kh_counter_t *hash, char *key);

// Clean up a counter
void kh_counter_destroy(kh_counter_t *hash);

// Increment a counter by one
size_t kh_counter_increment(kh_counter_t *hash, char *key);

// Decrement a counter by one
size_t kh_counter_sub(kh_counter_t *hash, char *key, int val);

// Increment a counter by a given amount
size_t kh_counter_add(kh_counter_t *hash, char *key, int val);

#endif
