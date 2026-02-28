#include "main.h"

#include <stdio.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "../fc_log.h"
#include "args.h"
#include "coverage.h"

int bamcoverage_main(int argc, char **argv)
{
    fc_log_config log_cfg = {
        .level = FC_LOG_INFO,
        .module = "BAMCOVER",
        .target = NULL,
        .color = FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    arguments_t args = bamcoverage_parse_arguments(argc, argv);

    htsThreadPool p = { NULL, 0 };
    p.pool = hts_tpool_init(args.threads);

    htsFile *in = sam_open(args.bam, "r");
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &p);

    uint32_t thresholds[] = { 1, 5, 10, 20 };
    size_t n_thresholds = sizeof(thresholds) / sizeof(thresholds[0]);

    cov_writer writer = init_coverage_writer(
        "bamstats-coverage", true, false,
        -1, -1, false,
        hdr, &p,
        args.beds, args.bed_names, args.n_beds,
        thresholds, n_thresholds,
        NULL, 0);
    bam1_t *rec = bam_init1();
    int ret;
    size_t records_read = 0;
    while ((ret = sam_read1(in, hdr, rec)) >= 0) {
        records_read++;
        coverage_process(writer, rec);
    }
    if (ret < -1) {
        FC_LOG_FATAL("Error reading BAM file (sam_read1 returned %d)", ret);
    }
    bam_destroy1(rec);
    destroy_coverage_writer(writer);

    bam_hdr_destroy(hdr);
    hts_close(in);

    hts_tpool_destroy(p.pool);

    FC_LOG_INFO("Processed %zu records from %s", records_read, args.bam);

    fc_log_destroy();
    return 0;
}
