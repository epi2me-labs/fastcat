#include "bamindex/main.h"

#include <stdio.h>
#include <string.h>

#include "fastlint/main.h"
#include "fc_catstats_args.h"
#include "fc_program.h"
#include "fc_program_kind.h"

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s <subcommand> [options]\n", prog);
    fprintf(stderr, "Subcommands:\n");
    fprintf(stderr, "  fastq   FASTQ concatenation/statistics/filtering mode (fastcat)\n");
    fprintf(stderr, "  bam     BAM stats/filtering mode (bamstats)\n");
    fprintf(stderr, "  lint    Standalone FASTQ low-complexity filter (fastlint)\n");
    fprintf(stderr, "  index   BAM chunk index tools (bamindex)\n");
    fprintf(stderr, "Try '%s <subcommand> --help' for subcommand options.\n", prog);
}

int main(int argc, char **argv)
{
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return (argc < 2) ? 1 : 0;
    }

    fc_program_kind program = FC_PROGRAM_FASTCAT;
    if (strcmp(argv[1], "fastq") == 0) {
        program = FC_PROGRAM_FASTCAT;
    } else if (strcmp(argv[1], "bam") == 0) {
        program = FC_PROGRAM_BAMSTATS;
    } else if (strcmp(argv[1], "lint") == 0 || strcmp(argv[1], "fastlint") == 0) {
        return fastlint_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "index") == 0 || strcmp(argv[1], "bamindex") == 0) {
        return bamindex_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "cover") == 0 || strcmp(argv[1], "bamcoverage") == 0) {
        fprintf(stderr, "Subcommand '%s' has been removed.\n", argv[1]);
        fprintf(stderr, "Use 'fastcat bam --coverage ...' instead.\n");
        return 1;
    } else {
        fprintf(stderr, "Unknown subcommand '%s'.\n", argv[1]);
        print_usage(argv[0]);
        return 1;
    }

    fc_catstats_args args = fc_parse_catstats_args(program, argc - 1, argv + 1);
    int status = fc_run_program(program, &args);
    fc_destroy_catstats_args(&args);
    return status;
}
