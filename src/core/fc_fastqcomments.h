#ifndef FASTCAT_FASTQCOMMENTS_H
#define FASTCAT_FASTQCOMMENTS_H

#include "htslib/kstring.h"
#include "htslib/sam.h"

int fastq_comment_to_bam(const kstring_t *comment, bam1_t *b);

#endif
