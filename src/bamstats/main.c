#include "fc_catstats_args.h"
#include "fc_program.h"
#include "fc_program_kind.h"

int main(int argc, char **argv) {
    fc_catstats_args args = fc_parse_catstats_args(FC_PROGRAM_BAMSTATS, argc, argv);
    int status = fc_run_program(FC_PROGRAM_BAMSTATS, &args);
    fc_destroy_catstats_args(&args);
    return status;
}
