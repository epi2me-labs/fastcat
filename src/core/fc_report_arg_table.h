#ifndef FASTCAT_FC_REPORT_ARG_TABLE_H
#define FASTCAT_FC_REPORT_ARG_TABLE_H

#include <stddef.h>

#include "cli/catstats/args.h"
#include "cli/catstats/program_kind.h"
#include "core/fc_report.h"

typedef struct {
    fc_report_config_entry *entries;
    size_t n_entries;
} fc_report_arg_table;

fc_report_arg_table fc_report_arg_table_from_catstats(
    fc_program_kind program,
    const fc_catstats_args *args,
    size_t file_count);
void fc_report_arg_table_destroy(fc_report_arg_table *table);

#endif
