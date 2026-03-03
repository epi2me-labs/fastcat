#ifndef FASTCAT_FC_PROGRAM_H
#define FASTCAT_FC_PROGRAM_H

#include "cli/catstats/args.h"
#include "cli/catstats/program_kind.h"

int fc_run_program(fc_program_kind program, const fc_catstats_args *args);

#endif
