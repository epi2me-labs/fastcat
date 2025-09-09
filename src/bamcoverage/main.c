#include <stdio.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "coverage.h"
#include "args.h"

int main(int argc, char** argv) {
    arguments_t args = parse_arguments(argc, argv);

    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(args.threads);
    
    htsFile* in = sam_open(args.bam, "r");
    bam_hdr_t* hdr = sam_hdr_read(in);
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &p);

    uint32_t thresholds[] = {1, 5, 10, 20};
    size_t n_thresholds = sizeof(thresholds) / sizeof(thresholds[0]);

    cov_writer writer = init_coverage_writer(
        "bamstats-coverage", true, false,
        -1, -1, false,
        hdr, &p,
        args.beds, args.bed_names, args.n_beds,
        thresholds, n_thresholds,
        NULL, 0);
    bam1_t* rec = bam_init1();
    while (sam_read1(in, hdr, rec) >= 0) {
        coverage_process(writer, rec);
    }
    bam_destroy1(rec);
    destroy_coverage_writer(writer);

    bam_hdr_destroy(hdr);
    hts_close(in);

    hts_tpool_destroy(p.pool);

    fprintf(stderr, "Processed %zu regions from %s\n", args.n_beds, args.bam);

    return 0;
}
