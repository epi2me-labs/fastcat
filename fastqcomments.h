#ifndef FASTCAT_FASTQCOMMENTS_H
#define FASTCAT_FASTQCOMMENTS_H

#include "kstring.h"

typedef struct {
    char* comment;
    char* runid;
    char* flow_cell_id;
    char* barcode;
    size_t ibarcode;
} read_meta;

read_meta parse_read_meta(kstring_t comment);

#endif
