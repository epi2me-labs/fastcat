#include "core/fc_alignment_summary.h"

#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_filter.h"

const char *const fc_flagstat_column_names[FC_FLAGSTAT_COUNT] = {
    "total",
    "primary",
    "secondary",
    "supplementary",
    "unmapped",
    "qcfail",
    "duplicate",
    "duplex",
    "duplex_forming"
};

struct fc_alignment_summary {
    bool region_mode;
    bool use_barcodes;
    size_t barcode_count;
    fc_alignment_counts *global;
    fc_alignment_counts **barcodes;
};

static fc_alignment_counts *alignment_counts_create(size_t n_refs)
{
    fc_alignment_counts *counts = fc_alloc(1, sizeof(*counts), "alignment counts");
    counts->n_refs = n_refs;
    counts->reference_counts = n_refs > 0
        ? fc_alloc(n_refs, sizeof(*counts->reference_counts), "alignment reference counts")
        : NULL;
    counts->unmapped_counts = fc_alloc(
        FC_FLAGSTAT_COUNT,
        sizeof(*counts->unmapped_counts),
        "alignment unmapped counts");
    for (size_t i = 0; i < n_refs; ++i) {
        counts->reference_counts[i] = fc_alloc(
            FC_FLAGSTAT_COUNT,
            sizeof(*counts->reference_counts[i]),
            "alignment reference count row");
    }
    return counts;
}

static void alignment_counts_destroy(fc_alignment_counts *counts)
{
    if (NULL == counts) {
        return;
    }
    for (size_t i = 0; i < counts->n_refs; ++i) {
        free(counts->reference_counts[i]);
    }
    free(counts->reference_counts);
    free(counts->unmapped_counts);
    free(counts);
}

fc_alignment_summary *fc_alignment_summary_create(
    size_t n_refs,
    bool region_mode,
    bool use_barcodes,
    size_t barcode_count)
{
    fc_alignment_summary *summary = fc_alloc(1, sizeof(*summary), "alignment summary");
    summary->region_mode = region_mode;
    summary->use_barcodes = use_barcodes;
    summary->barcode_count = use_barcodes && barcode_count > 0 ? barcode_count : 1;
    size_t stored_refs = region_mode ? 1 : n_refs;
    summary->global = alignment_counts_create(stored_refs);
    if (use_barcodes) {
        summary->barcodes = fc_alloc(
            summary->barcode_count,
            sizeof(*summary->barcodes),
            "alignment barcode counts");
    }
    return summary;
}

void fc_alignment_summary_destroy(fc_alignment_summary *summary)
{
    if (NULL == summary) {
        return;
    }
    alignment_counts_destroy(summary->global);
    for (size_t i = 0; i < summary->barcode_count && NULL != summary->barcodes; ++i) {
        alignment_counts_destroy(summary->barcodes[i]);
    }
    free(summary->barcodes);
    free(summary);
}

static size_t *alignment_counts_row(
    const fc_alignment_summary *summary,
    fc_alignment_counts *counts,
    const bam1_t *bam)
{
    if (fc_is_unmapped(bam)) {
        return counts->unmapped_counts;
    }
    if (0 == counts->n_refs || NULL == counts->reference_counts) {
        return NULL;
    }
    size_t tid = (size_t)bam->core.tid;
    if (!summary->region_mode && tid >= counts->n_refs) {
        return NULL;
    }
    return summary->region_mode ? counts->reference_counts[0] : counts->reference_counts[tid];
}

static void alignment_counts_add(size_t *counts, const bam1_t *bam, int duplex_code)
{
    counts[FC_FLAGSTAT_TOTAL] += 1;
    counts[FC_FLAGSTAT_PRIMARY] += fc_is_primary(bam);
    counts[FC_FLAGSTAT_SECONDARY] += ((bam->core.flag & BAM_FSECONDARY) != 0);
    counts[FC_FLAGSTAT_SUPPLEMENTARY] += ((bam->core.flag & BAM_FSUPPLEMENTARY) != 0);
    counts[FC_FLAGSTAT_UNMAPPED] += ((bam->core.flag & BAM_FUNMAP) != 0);
    counts[FC_FLAGSTAT_QCFAIL] += ((bam->core.flag & BAM_FQCFAIL) != 0);
    counts[FC_FLAGSTAT_DUPLICATE] += ((bam->core.flag & BAM_FDUP) != 0);
    counts[FC_FLAGSTAT_DUPLEX] += (duplex_code == 1);
    counts[FC_FLAGSTAT_DUPLEX_FORMING] += (duplex_code == -1);
}

static bool alignment_counts_add_read(
    const fc_alignment_summary *summary,
    fc_alignment_counts *counts,
    const bam1_t *bam,
    int duplex_code)
{
    size_t *row = alignment_counts_row(summary, counts, bam);
    if (NULL == row) {
        return false;
    }
    alignment_counts_add(row, bam, duplex_code);
    counts->seen = true;
    return true;
}

void fc_alignment_summary_add_read(
    fc_alignment_summary *summary,
    size_t barcode,
    const bam1_t *bam,
    int duplex_code)
{
    if (NULL == summary || NULL == bam) {
        return;
    }
    if (!alignment_counts_add_read(summary, summary->global, bam, duplex_code)) {
        return;
    }
    if (!summary->use_barcodes || barcode >= summary->barcode_count) {
        return;
    }
    if (NULL == summary->barcodes[barcode]) {
        summary->barcodes[barcode] = alignment_counts_create(summary->global->n_refs);
    }
    alignment_counts_add_read(summary, summary->barcodes[barcode], bam, duplex_code);
}

const fc_alignment_counts *fc_alignment_summary_global(
    const fc_alignment_summary *summary)
{
    return NULL != summary ? summary->global : NULL;
}

const fc_alignment_counts *fc_alignment_summary_barcode(
    const fc_alignment_summary *summary,
    size_t barcode)
{
    if (NULL == summary
        || !summary->use_barcodes
        || barcode >= summary->barcode_count
        || NULL == summary->barcodes) {
        return NULL;
    }
    return summary->barcodes[barcode];
}

const fc_alignment_counts *fc_alignment_summary_for_output(
    const fc_alignment_summary *summary,
    size_t barcode)
{
    if (NULL == summary) {
        return NULL;
    }
    return summary->use_barcodes
        ? fc_alignment_summary_barcode(summary, barcode)
        : fc_alignment_summary_global(summary);
}

void fc_alignment_counts_materialize(
    const fc_alignment_counts *counts,
    fc_alignment_totals *totals)
{
    if (NULL == totals) {
        return;
    }
    memset(totals, 0, sizeof(*totals));
    if (NULL == counts) {
        return;
    }
    for (size_t ref = 0; ref < counts->n_refs; ++ref) {
        for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
            totals->record_counts[col] += counts->reference_counts[ref][col];
        }
        totals->mapped_record_count += counts->reference_counts[ref][FC_FLAGSTAT_TOTAL];
        totals->primary_mapped_record_count += counts->reference_counts[ref][FC_FLAGSTAT_PRIMARY];
    }
    for (size_t col = 0; col < FC_FLAGSTAT_COUNT; ++col) {
        totals->record_counts[col] += counts->unmapped_counts[col];
    }
    totals->read_count = totals->primary_mapped_record_count
        + totals->record_counts[FC_FLAGSTAT_UNMAPPED];
    totals->mapping_rate = totals->read_count > 0
        ? (double)totals->primary_mapped_record_count / (double)totals->read_count
        : 0.0;
}
