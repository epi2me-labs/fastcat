// Wrap khash to make it more consise to use

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include "kh_counter.h"

/* Implementation of a counter of strings (increasing only)
 *
 * kh_counter_t *counter = kh_counter_init();
 * kh_counter_increment(counter, "one");
 * kh_counter_increment(counter, "two");
 * kh_counter_increment(counter, "two");
 * kh_counter_add(counter, "three", 2);
 * kh_counter_increment(counter, "three");
 * kh_counter_destroy(h);
 *
 */


kh_counter_t *kh_counter_init(void) {
    kh_counter_t *h = kh_init(KH_COUNTER);
    return h;
}

int kh_counter_val(kh_counter_t *hash, char *key) {
    khiter_t k = kh_get(KH_COUNTER, hash, key);
    int val = k != kh_end(hash) ? kh_val(hash, k) : 0;
    return val;
}

size_t kh_counter_add(kh_counter_t *hash, char *key, int val) {
    if (key == NULL) {return -1;}
    // note: key is copied so no need for caller to hold on to it
    int ret;
    khiter_t k = kh_put(KH_COUNTER, hash, key, &ret);
    if (ret == 1) { // new key
        kh_key(hash, k) = strdup(key);
        kh_value(hash, k) = val;
    } else if (ret == 0) {  // exists
        // get value and add
        int cur = kh_val(hash, k);
        kh_value(hash, k) = cur + val;
    } else {
        // shouldnt get here - previously deleted key
    }
    return ret;
}

size_t kh_counter_sub(kh_counter_t *hash, char *key, int val) {
    if (key == NULL) {return -1;}
    // note: key is copied so no need for caller to hold on to it
    int ret;
    khiter_t k = kh_put(KH_COUNTER, hash, key, &ret);
    if (ret == 1) { // new key
        kh_key(hash, k) = strdup(key);
        kh_value(hash, k) = -val;
    } else if (ret == 0) {  // exists
        // get value and add
        int cur = kh_val(hash, k);
        kh_value(hash, k) = cur - val;
    } else {
        // shouldnt get here - previously deleted key
    }
    return ret;
}

size_t kh_counter_increment(kh_counter_t *hash, char *key) {
    return kh_counter_add(hash, key, 1);
}

void kh_counter_destroy(kh_counter_t *hash) {
    for (khiter_t k = 0; k < kh_end(hash); k++){
        if (kh_exist(hash, k)) {
            free((char*) kh_key(hash, k));
        }
    }
    kh_destroy(KH_COUNTER, hash);
}
