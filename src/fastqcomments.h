#ifndef FASTCAT_FASTQCOMMENTS_H
#define FASTCAT_FASTQCOMMENTS_H

#include "common.h"
#include "htslib/kstring.h"

typedef struct {
    char* comment;
    char* rg;
    readgroup* rg_info;
    char* runid;
    char* basecaller;
    char* flow_cell_id;
    char* barcode;
    size_t ibarcode;
    char* barcode_alias;
    char* start_time;
    size_t read_number;
    size_t channel;
    kstring_t* rest;
    kstring_t* tags_str;
} _read_meta;

typedef _read_meta* read_meta;


// constructor
read_meta create_read_meta(const kstring_t* comment);

// destructor
void destroy_read_meta(read_meta meta);

// parser
read_meta parse_read_meta(kstring_t comment);

#endif
