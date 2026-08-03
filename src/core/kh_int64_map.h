#ifndef FASTCAT_KH_INT64_MAP_H
#define FASTCAT_KH_INT64_MAP_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "htslib/khash.h"

KHASH_MAP_INIT_INT64(KH_INT64_MAP, void *)

typedef khash_t(KH_INT64_MAP) kh_int64_map;
typedef void (*kh_int64_map_visit_fn)(int64_t key, void *value, void *context);

kh_int64_map *kh_int64_map_create(void);
void kh_int64_map_destroy(kh_int64_map *map, kh_int64_map_visit_fn visit, void *context);

void *kh_int64_map_get(const kh_int64_map *map, int64_t key);
bool kh_int64_map_insert(kh_int64_map *map, int64_t key, void *value);
size_t kh_int64_map_size(const kh_int64_map *map);
void kh_int64_map_foreach(
    const kh_int64_map *map, kh_int64_map_visit_fn visit, void *context);

#endif
