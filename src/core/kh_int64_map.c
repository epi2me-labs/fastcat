#include "core/kh_int64_map.h"

#include <stdlib.h>

kh_int64_map *kh_int64_map_create(void)
{
    kh_int64_map *map = kh_init(KH_INT64_MAP);
    return map;
}

void kh_int64_map_foreach(const kh_int64_map *map, kh_int64_map_visit_fn visit, void *context)
{
    if (!map || !visit) {
        return;
    }
    for (khiter_t i = kh_begin(map); i != kh_end(map); ++i) {
        if (kh_exist(map, i)) {
            visit((int64_t)kh_key(map, i), kh_val(map, i), context);
        }
    }
}

void kh_int64_map_destroy(kh_int64_map *map, kh_int64_map_visit_fn visit, void *context)
{
    if (!map) {
        return;
    }
    kh_int64_map_foreach(map, visit, context);
    kh_destroy(KH_INT64_MAP, map);
}

void *kh_int64_map_get(const kh_int64_map *map, int64_t key)
{
    if (!map) {
        return NULL;
    }
    khiter_t i = kh_get(KH_INT64_MAP, map, (khint64_t)key);
    return i == kh_end(map) ? NULL : kh_val(map, i);
}

bool kh_int64_map_insert(kh_int64_map *map, int64_t key, void *value)
{
    if (!map) {
        return false;
    }
    int ret = 0;
    khiter_t i = kh_put(KH_INT64_MAP, map, (khint64_t)key, &ret);
    if (ret < 0 || ret == 0) {
        return false;
    }
    kh_val(map, i) = value;
    return true;
}

size_t kh_int64_map_size(const kh_int64_map *map)
{
    return map ? kh_size(map) : 0;
}
