#ifndef _BAMSTATS_STATS_H
#define _BAMSTATS_STATS_H

#include <stdbool.h>
#include "htslib/sam.h"

#include "args.h"
#include "../bamcoverage/coverage.h"
#include "../fc_stats_engine.h"
#include "../kh_counter.h"
#include "../common.h"

// struct for flagstat counts
typedef struct flag_stats {
    size_t n_refs;
    size_t** counts;
    size_t* unmapped;
} flag_stats;

/** Create flagstat counts struct for a BAM file.
 *
 *  @param n_refs number of reference sequences.
 *  @param store_unmapped whether to count unmapped reads.
 *
 */
flag_stats* create_flag_stats(size_t n_refs, bool store_unmapped);

/** Clean up flagstat counts.
 *
 *  @param stats flagstat counts structure to clean.
 *
 */
void destroy_flag_stats(flag_stats* stats);


#endif
