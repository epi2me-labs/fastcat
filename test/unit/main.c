#include "unity.h"

void common_suite(void);
void stats_suite(void);
void fastqcomments_suite(void);
void filter_suite(void);
void rg_parse_suite(void);
void regiter_suite(void);
void quality_suite(void);
void readstats_suite(void);
void bam_utils_suite(void);

void setUp(void) {}
void tearDown(void) {}

int main(void) {
    UNITY_BEGIN();
    common_suite();
    stats_suite();
    fastqcomments_suite();
    filter_suite();
    rg_parse_suite();
    regiter_suite();
    quality_suite();
    readstats_suite();
    bam_utils_suite();
    return UNITY_END();
}
