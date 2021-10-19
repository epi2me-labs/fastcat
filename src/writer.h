#ifndef FASTCAT_WRITER_H
#define FASTCAT_WRITER_H

#include <zlib.h>

#include "kseq.h"
#include "fastqcomments.h"

KSEQ_INIT(gzFile, gzread)

// barcode 0 is reserved for "unclassified"
#define MAX_BARCODES 1025

typedef struct {
    char* path;
    char* output;
    gzFile* handles;
    size_t* nreads;
    char* sample;
    size_t reheader;
    FILE* perread;
    FILE* perfile;
} _writer;

typedef _writer* writer;

char* strip_path(char* input);

writer initialize_writer(char* output_dir, char* perread, char* perfile, char* sample, size_t reheader);

void destroy_writer(writer writer);

void write_read(writer writer, kseq_t* seq, read_meta meta, float mean_q, char* fname);

#endif
