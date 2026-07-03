#include "unity.h"
#include "core/fc_log.h"

void common_suite(void);
void stats_suite(void);
void stats_engine_suite(void);
void fastqcomments_suite(void);
void filter_suite(void);
void rg_parse_suite(void);
void regiter_suite(void);
void quality_suite(void);
void bam_utils_suite(void);
void qualcalcs_suite(void);
void iter_suite(void);
void iter_collect_suite(void);
void bamcoverage_suite(void);
void quantile_suite(void);
void kll_suite(void);
void length_qscore_2d_suite(void);
void time_bin_suite(void);
void quantile_stress_suite(void);

void setUp(void) {}
void tearDown(void) {}

int main(void) {
    fc_log_config log_cfg = {
        .level = FC_LOG_DEBUG,
        .module = "UNIT_TEST",
        .target = NULL,
        .color = FC_LOG_COLOR_OFF
    };
    fc_log_init(&log_cfg);

    UNITY_BEGIN();
    common_suite();
    stats_suite();
    stats_engine_suite();
    fastqcomments_suite();
    filter_suite();
    rg_parse_suite();
    regiter_suite();
    quality_suite();
    bam_utils_suite();
    qualcalcs_suite();
    iter_suite();
    iter_collect_suite();
    bamcoverage_suite();
    quantile_suite();
    kll_suite();
    length_qscore_2d_suite();
    time_bin_suite();
    quantile_stress_suite();
    int rc = UNITY_END();

    fc_log_destroy();
    return rc;
}
