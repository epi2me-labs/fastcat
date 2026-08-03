#ifndef FASTCAT_FC_ALIGNMENT_SUMMARY_H
#define FASTCAT_FC_ALIGNMENT_SUMMARY_H

#include <stdbool.h>
#include <stddef.h>

#include "htslib/sam.h"

typedef enum {
    FC_FLAGSTAT_TOTAL,
    FC_FLAGSTAT_PRIMARY,
    FC_FLAGSTAT_SECONDARY,
    FC_FLAGSTAT_SUPPLEMENTARY,
    FC_FLAGSTAT_UNMAPPED,
    FC_FLAGSTAT_QCFAIL,
    FC_FLAGSTAT_DUPLICATE,
    FC_FLAGSTAT_DUPLEX,
    FC_FLAGSTAT_DUPLEX_FORMING,
    FC_FLAGSTAT_COUNT
} fc_flagstat_column;

extern const char *const fc_flagstat_column_names[FC_FLAGSTAT_COUNT];

typedef struct fc_alignment_summary fc_alignment_summary;

typedef struct {
    bool seen;
    size_t n_refs;
    // Per-reference mapped-record counts followed by a separate unmapped row.
    size_t **reference_counts;
    size_t *unmapped_counts;
    size_t canonical_unmapped_read_count;
} fc_alignment_counts;

typedef struct {
    size_t record_counts[FC_FLAGSTAT_COUNT];
    size_t mapped_record_count;
    size_t primary_mapped_record_count;
    // Read-level totals use canonical mapped and unmapped records only.
    size_t mapped_read_count;
    size_t read_count;
    double mapping_rate;
} fc_alignment_totals;

fc_alignment_summary *fc_alignment_summary_create(
    size_t n_refs,
    bool region_mode,
    bool use_barcodes,
    size_t barcode_count);
void fc_alignment_summary_destroy(fc_alignment_summary *summary);

void fc_alignment_summary_add_read(
    fc_alignment_summary *summary,
    size_t barcode,
    const bam1_t *bam,
    int duplex_code);

const fc_alignment_counts *fc_alignment_summary_global(
    const fc_alignment_summary *summary);
const fc_alignment_counts *fc_alignment_summary_barcode(
    const fc_alignment_summary *summary,
    size_t barcode);
const fc_alignment_counts *fc_alignment_summary_for_output(
    const fc_alignment_summary *summary,
    size_t barcode);

void fc_alignment_counts_materialize(
    const fc_alignment_counts *counts,
    fc_alignment_totals *totals);

#endif
