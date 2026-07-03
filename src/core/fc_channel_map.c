#include "core/fc_channel_map.h"

void fc_channel_map_project(const fc_channel_map *map, double raw_x, double raw_y, double *plot_x, double *plot_y)
{
    if (NULL == map || NULL == plot_x || NULL == plot_y) {
        return;
    }
    if (map->swap_axes) {
        *plot_x = raw_y;
        *plot_y = -raw_x;
        return;
    }
    *plot_x = raw_x;
    *plot_y = -raw_y;
}

fc_channel_map_device fc_channel_map_device_for_channel(size_t highest_channel)
{
    if (highest_channel <= 128) {
        return FC_CHANNEL_MAP_FLONGLE;
    }
    if (highest_channel <= 512) {
        return FC_CHANNEL_MAP_MINION;
    }
    return FC_CHANNEL_MAP_PROMETHION;
}

const fc_channel_map *fc_channel_map_get(fc_channel_map_device device)
{
    static fc_channel_map maps[] = {
        { .rows = fc_channel_map_flongle_rows,
            .row_count = 0,
            .x_min = 0.0,
            .x_max = 0.0,
            .y_min = 0.0,
            .y_max = 0.0 },
        { .rows = fc_channel_map_minion_rows,
            .row_count = 0,
            .x_min = 0.0,
            .x_max = 0.0,
            .y_min = 0.0,
            .y_max = 0.0 },
        { .rows = fc_channel_map_promethion_rows,
            .row_count = 0,
            .x_min = 0.0,
            .x_max = 0.0,
            .y_min = 0.0,
            .y_max = 0.0 }
    };
    static int initialized = 0;

    if (!initialized) {
        maps[FC_CHANNEL_MAP_FLONGLE].row_count = fc_channel_map_flongle_row_count;
        maps[FC_CHANNEL_MAP_FLONGLE].swap_axes = true;
        maps[FC_CHANNEL_MAP_FLONGLE].x_min = (double)fc_channel_map_flongle_y_min;
        maps[FC_CHANNEL_MAP_FLONGLE].x_max = (double)fc_channel_map_flongle_y_max;
        maps[FC_CHANNEL_MAP_FLONGLE].y_min = -(double)fc_channel_map_flongle_x_max;
        maps[FC_CHANNEL_MAP_FLONGLE].y_max = -(double)fc_channel_map_flongle_x_min;

        maps[FC_CHANNEL_MAP_MINION].row_count = fc_channel_map_minion_row_count;
        maps[FC_CHANNEL_MAP_MINION].swap_axes = true;
        maps[FC_CHANNEL_MAP_MINION].x_min = (double)fc_channel_map_minion_y_min;
        maps[FC_CHANNEL_MAP_MINION].x_max = (double)fc_channel_map_minion_y_max;
        maps[FC_CHANNEL_MAP_MINION].y_min = -(double)fc_channel_map_minion_x_max;
        maps[FC_CHANNEL_MAP_MINION].y_max = -(double)fc_channel_map_minion_x_min;

        maps[FC_CHANNEL_MAP_PROMETHION].row_count = fc_channel_map_promethion_row_count;
        maps[FC_CHANNEL_MAP_PROMETHION].swap_axes = true;
        maps[FC_CHANNEL_MAP_PROMETHION].x_min = (double)fc_channel_map_promethion_y_min;
        maps[FC_CHANNEL_MAP_PROMETHION].x_max = (double)fc_channel_map_promethion_y_max;
        maps[FC_CHANNEL_MAP_PROMETHION].y_min = -(double)fc_channel_map_promethion_x_max;
        maps[FC_CHANNEL_MAP_PROMETHION].y_max = -(double)fc_channel_map_promethion_x_min;
        initialized = 1;
    }

    if (device > FC_CHANNEL_MAP_PROMETHION) {
        return NULL;
    }
    return &maps[device];
}
