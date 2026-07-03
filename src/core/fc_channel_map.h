#ifndef FASTCAT_FC_CHANNEL_MAP_H
#define FASTCAT_FC_CHANNEL_MAP_H

#include <stdbool.h>
#include <stddef.h>

#include "core/fc_channel_map_data.h"

typedef enum {
    FC_CHANNEL_MAP_FLONGLE = 0,
    FC_CHANNEL_MAP_MINION,
    FC_CHANNEL_MAP_PROMETHION
} fc_channel_map_device;

typedef struct {
    const fc_channel_map_channel_row *rows;
    size_t row_count;
    bool swap_axes;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
} fc_channel_map;

fc_channel_map_device fc_channel_map_device_for_channel(size_t highest_channel);
const fc_channel_map *fc_channel_map_get(fc_channel_map_device device);
void fc_channel_map_project(const fc_channel_map *map, double raw_x, double raw_y, double *plot_x, double *plot_y);

#endif
