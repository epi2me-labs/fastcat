#ifndef FASTCAT_FC_CHANNEL_MAP_DATA_H
#define FASTCAT_FC_CHANNEL_MAP_DATA_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint16_t channel;
    int16_t x[4];
    int16_t y[4];
} fc_channel_map_channel_row;

extern const fc_channel_map_channel_row fc_channel_map_flongle_rows[];
extern const size_t fc_channel_map_flongle_row_count;
extern const int fc_channel_map_flongle_x_min;
extern const int fc_channel_map_flongle_x_max;
extern const int fc_channel_map_flongle_y_min;
extern const int fc_channel_map_flongle_y_max;

extern const fc_channel_map_channel_row fc_channel_map_minion_rows[];
extern const size_t fc_channel_map_minion_row_count;
extern const int fc_channel_map_minion_x_min;
extern const int fc_channel_map_minion_x_max;
extern const int fc_channel_map_minion_y_min;
extern const int fc_channel_map_minion_y_max;

extern const fc_channel_map_channel_row fc_channel_map_promethion_rows[];
extern const size_t fc_channel_map_promethion_row_count;
extern const int fc_channel_map_promethion_x_min;
extern const int fc_channel_map_promethion_x_max;
extern const int fc_channel_map_promethion_y_min;
extern const int fc_channel_map_promethion_y_max;

#endif
