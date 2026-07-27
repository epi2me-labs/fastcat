#ifndef FASTCAT_FC_READ_H
#define FASTCAT_FC_READ_H

#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <stdbool.h>

#include "core/fc_bam_rg_cache.h"
#include "core/fc_common.h"
#include "core/fc_stats.h"
#include "core/kseq_gz.h"

typedef struct read_meta {
    char *orig_comment;
    char *rg; // RG:Z
    readgroup *rg_info;
    char *runid; // RG:Z (RD:Z legacy)
    char *basecaller;
    char *_basecaller; // owned storage for basecaller when needed
    char *modcallers;
    char *_modcallers;
    char *flow_cell_id; // PU:Z
    char *sample_id; // LB:Z
    char *barcode; // SM:Z
    size_t ibarcode;
    char *barcode_alias; // al:Z
    char *start_time; // st:Z
    int edit_distance; // NM:i
    char *parent_read; // pi:Z (parent_read_id)
    int poly_tail_length; // pt:i
    int duplex_code; // dx:i
    float mean_q; // qs:f
    double duration_s; // du:f
    size_t read_number; // rn:i
    size_t channel; // ch:i
    size_t mux; // mx:i
    kstring_t *tags_str; // lazily materialized string form of aux tags
    bool tags_ready;
} read_meta;

typedef struct {
    bam1_t *bam;
    read_meta *meta;
} fc_read;

// Allocate and initialize a heap-owned fc_read and its BAM record storage.
fc_read *fc_read_init(void);

// Destroy a caller-owned fc_read and free all owned allocations.
void fc_read_destroy(fc_read *read);

// Convert a FASTQ kseq record into an unaligned BAM record, including FASTQ comment tags.
int kseq_to_bam(bam1_t *b, kseq_t *seq);

// Reset a reusable fc_read in-place from the current BAM record.
// Expects mean quality tag policy to have been applied by the caller.
void fc_read_reset_from_bam(fc_read *read, const bam_hdr_t *hdr, const fc_bam_rg_cache *rg_cache, const kstring_t *orig_comment);

// Return the demultiplexing barcode index without constructing full read metadata.
// This is used while accounting for reads rejected before normal metadata parsing.
size_t fc_read_barcode_index_from_bam(
    const bam1_t *bam,
    const bam_hdr_t *hdr,
    const fc_bam_rg_cache *rg_cache);

// Return SAM-style auxiliary tags text, lazily materialized and cached on first use.
const kstring_t *fc_read_get_tags_str(const fc_read *read);

#endif
