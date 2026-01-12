#ifndef FASTCAT_BAM_UTILS_H
#define FASTCAT_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "../fastqcomments.h"

#ifndef KSEQ_DECLARED
#include <zlib.h>
KSEQ_DECLARE(gzFile)
#define KSEQ_DECLARED
#endif

// Parse a tab-delimited SAM-style aux string and attach tags to a bam1_t
int parse_and_set_aux_tags(bam1_t *b, const char *aux_str);

// Convert a FASTQ kseq record plus parsed metadata into an unmapped bam1_t
// Caller owns the returned bam1_t and must free with bam_destroy1.
bam1_t* kseq_to_bam(kseq_t* seq, read_meta meta);

#endif
