#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stats.h"
#include "common.h"


read_stats* create_length_stats(void) {
    read_stats* stats = (read_stats*) xalloc(1, sizeof(read_stats), "length_stats");

    bin_groups* bins = (bin_groups*) xalloc(1, sizeof(bin_groups), "bin_groups");
    stats->buckets = bins;
    bins->n = 1;
    bins->groups = (size_t*) xalloc(1, 3*bins->n*sizeof(size_t), "groups");
    size_t* grps = bins->groups;
    // - end exclusive upper edges
    // - 0 is first lower edge
    // - final >x bucket
    //grps[0] =   50000; grps[1] =     1;
    //grps[3] =  100000; grps[4] =    10;
    //grps[6] = 1000000; grps[7] = 1000;
    // just do one massive bucket up to 10M - 76Mbytes
    grps[0] = 10000000; grps[1] = 1;


    // count the total number of bins across all groups
    stats->n = 0;
    size_t lower = 0;
    for (size_t i=0; i<bins->n; i++) {
        size_t upper = bins->groups[3*i];
        size_t step = bins->groups[3*i + 1];
        size_t nbins = (upper - lower) / step;
        bins->groups[3*i + 2] = nbins;
        stats->n += nbins;
        lower = upper;
    }
    stats->n++;

    // fill in all the edges
    stats->width = 0;
    stats->edges = xalloc(stats->n, sizeof(size_t), "edges");
    stats->counts = xalloc(stats->n, sizeof(size_t), "counts");
    size_t i=0;
    lower = 0;
    size_t upper = 0;
    for (size_t b=0; b<bins->n; b++) {
        upper = bins->groups[3*b];
        size_t step = bins->groups[3*b + 1];
        for (size_t j=lower; j<upper; j+=step, i++) {
            stats->edges[i] = j;
        }
        lower = upper;
    }
    stats->edges[i] = upper;
    return stats;
}

void destroy_length_stats(read_stats* stats) {
    free(stats->buckets->groups);
    free(stats->buckets);
    free(stats->edges);
    free(stats->counts);
    free(stats);
}

void add_length_count(read_stats* stats, size_t x) {
    size_t lower = 0;
    size_t cum_bin = 0;
    bool done = false;
    for (size_t i=0; i<stats->buckets->n; i++) {
        size_t upper = stats->buckets->groups[3*i];
        if (x < upper) {
            size_t step = stats->buckets->groups[3*i + 1];
            stats->counts[cum_bin + (x - lower) / step]++;
            done = true;
            break;
        }
        lower = upper;
        cum_bin += stats->buckets->groups[3*i + 2];
    }
    if (!done) {
        stats->counts[cum_bin]++;
    }
}


read_stats* create_qual_stats(float width) {
    read_stats* stats = (read_stats*) xalloc(1, sizeof(read_stats), "quality stats");
    stats->width = width;
    // this fixes the range to [0, 100], good for both QUAL and %age acc
    stats->n = (size_t) (100.0 / stats->width) + 1;
    stats->counts = xalloc(stats->n, sizeof(size_t), "counts");
    return stats;
}

void destroy_qual_stats(read_stats* stats) {
    free(stats->counts);
    free(stats);
}

void add_qual_count(read_stats* stats, float q) {
    q = fmin(q, 100.0);
    stats->counts[(int) (q / stats->width)]++;
}

void print_stats(read_stats* stats, bool zeroes, bool tsv, FILE* fp) {
    if (fp == NULL) {
        fp = stderr;
    }
    if (stats->width == 0) {
        for (size_t i=0; i<stats->n; i++) {
            if (stats->counts[i] == 0 && !zeroes) continue;
            if (tsv) {
                fprintf(fp, "%zu\t%zu\t%zu\n", stats->edges[i], stats->edges[i+1], stats->counts[i]);
            }
            else {
                fprintf(fp, "[%zu, %zu)\t%zu\n", stats->edges[i], stats->edges[i+1], stats->counts[i]);
            }
        }
    }
    else {
        for (size_t i=0; i<stats->n; i++) {
            if (stats->counts[i] == 0 && !zeroes) continue;
            if (tsv) {
                fprintf(fp, "%.2f\t%.2f\t%zu\n", (float) i * stats->width, (float) (i+1) * stats->width, stats->counts[i]);
            }
            else {
                fprintf(fp, "[%.2f, %.2f)\t%zu\n", (float) i * stats->width, (float) (i+1) * stats->width, stats->counts[i]);
            }
        }
    }
}

//int main(int argc, char **argv) {
//    read_stats* stats = create_length_stats();
//    // bins every 1
//    add_length_count(stats, 1);
//    add_length_count(stats, 4);
//    add_length_count(stats, 950);
//    add_length_count(stats, 998);
//    add_length_count(stats, 999);
//
//    //// changing to bins every 10
//    add_length_count(stats, 1000);
//    add_length_count(stats, 1001);
//    add_length_count(stats, 1009);
//    add_length_count(stats, 1010);
//    add_length_count(stats, 1045);
//    add_length_count(stats, 1050);
//
//    // changing to bins every 100
//    add_length_count(stats, 9845);
//    add_length_count(stats, 9900);
//    add_length_count(stats, 9901);
//    add_length_count(stats, 9909);
//    add_length_count(stats, 9910);
//    add_length_count(stats, 9999);
//    add_length_count(stats, 10000);
//    add_length_count(stats, 10001);
//    add_length_count(stats, 10010);
//    add_length_count(stats, 10100);
//    add_length_count(stats, 10150);
//    add_length_count(stats, 10199);
//    add_length_count(stats, 10200);
//    add_length_count(stats, 10210);
//
//    // changing to bins every 1000
//    add_length_count(stats, 99100); 
//    add_length_count(stats, 99150); 
//    add_length_count(stats, 99500); 
//    add_length_count(stats, 99999); 
//    add_length_count(stats, 100000); 
//    add_length_count(stats, 100001); 
//    add_length_count(stats, 100999); 
//    add_length_count(stats, 101000); 
//    add_length_count(stats, 102000); 
//
//    add_length_count(stats,  999999);
//    add_length_count(stats, 1000000);
//    add_length_count(stats, 2000000);
//
//    print_stats(stats, false, false, stderr);
//    destroy_length_stats(stats);
//
//    float width = 0.02;
//    read_stats* qstats = create_q_stats(width);
//    
//    add_q_count(qstats, 9.89);
//    add_q_count(qstats, 9.89);
//    add_q_count(qstats, 99.5);
//    add_q_count(qstats, 99.9);
//    add_q_count(qstats, 100);
//    add_q_count(qstats, 101);
//    add_q_count(qstats, 101);
//    add_q_count(qstats, 101);
//    print_stats(qstats, false, false, stderr);
//    destroy_q_stats(qstats);
//
//}
