#ifndef FASTCAT_STATS_H
#define FASTCAT_STATS_H

#include "stdbool.h"
// use only 1 non-zero decimal place in these -- see _leading_decimals()
#define QUAL_HIST_WIDTH 0.02   // QUAL is log10
#define ACC_HIST_WIDTH 0.0001  // ACC is linear %age, Q60
#define COV_HIST_WIDTH 0.01    // COV is linear %age

typedef struct {
    size_t n;
    size_t* groups;  // 3 items: upper, step, number 
} bin_groups;


typedef struct {
    size_t n;
    float width;  // for fixed width
    size_t* edges;
    size_t* counts;
    bin_groups* buckets;
} read_stats;


read_stats* create_length_stats(void);
void destroy_length_stats(read_stats* stats);
void add_length_count(read_stats* stats, size_t);

// uses read_stats with a fixed (and rescaled grid)
read_stats* create_qual_stats(float width);
void destroy_qual_stats(read_stats* stats);
void add_qual_count(read_stats* stats, float q);

void print_stats(read_stats* stats, bool zeroes, bool tsv, FILE* fp);

size_t _leading_decimals(float num);
#endif
