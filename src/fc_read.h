#ifndef FASTCAT_FC_READ_H
#define FASTCAT_FC_READ_H

#include <stdbool.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include "common.h"
#include "stats.h"

#ifndef KSEQ_DECLARED
#include <zlib.h>
KSEQ_DECLARE(gzFile)
#define KSEQ_DECLARED
#endif

typedef struct {
    char* orig_comment;
    char* rg;             // RG:Z
    readgroup* rg_info;
    char* runid;          // RG:Z (RD:Z legacy)
    char* basecaller;
    char* flow_cell_id;   // PU:Z
    char* barcode;        // SM:Z
    size_t ibarcode;
    char* barcode_alias;  // al:Z
    char* start_time;     // st:Z
    int edit_distance;    // NM:i
    char* parent_read;    // pi:Z (parent_read_id)
    int poly_tail_length; // pt:i
    int duplex_code;      // dx:i
    float mean_q;         // qs:f
    size_t read_number;   // rn:i
    size_t channel;       // ch:i
    kstring_t* tags_str;
} _read_meta;

typedef _read_meta* read_meta;

typedef struct {
    bam1_t *bam;
    read_meta meta;
    bool owns_bam;
} fc_read;

// Resolve mean quality, preferring qs tag unless force_recalc_qual is true.
static inline float resolve_mean_quality(const bam1_t *b, const read_meta meta, bool force_recalc_qual) {
    float mean_quality = meta ? meta->mean_q : -1.0f;
    if (mean_quality == -1.0f || force_recalc_qual) {
        mean_quality = mean_qual_from_bam_naive(bam_get_qual(b), b->core.l_qseq);
    }
    return mean_quality;
}

// Free a read_meta and any associated allocations.
void destroy_read_meta(read_meta meta);

// Construct from an existing BAM record; does not take ownership of bam.
fc_read* fc_read_from_bam(bam1_t* bam, const bam_hdr_t* hdr, bool force_recalc_qual);

// Construct from FASTQ (kseq); owns the BAM.
fc_read* fc_read_from_kseq(kseq_t* seq);

// Free metadata/tags, and BAM if owned.
void fc_read_destroy(fc_read* r);

// Construct read_meta from BAM aux tags, optionally preserving the original FASTQ comment.
// Caller owns the returned read_meta and must destroy it.
read_meta create_read_meta(const bam1_t *b, const bam_hdr_t *hdr, bool require_nm,
                           const kstring_t *orig_comment);

#endif
