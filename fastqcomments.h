#ifndef FASTCAT_FASTQCOMMENTS_H
#define FASTCAT_FASTQCOMMENTS_H

#include "kstring.h"

typedef struct {
    char* comment;
    char* runid;
    char* flow_cell_id;
    char* barcode;
    size_t ibarcode;
    char* barcode_alias;
    char* start_time;
    size_t read_number;
    size_t channel;
} _read_meta;

typedef _read_meta* read_meta;

read_meta parse_read_meta(kstring_t comment);

void destroy_read_meta(read_meta meta);

#endif
