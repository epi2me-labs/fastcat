#ifndef FASTCAT_WRITER_H
#define FASTCAT_WRITER_H

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_BARCODES 8

typedef struct {
    char* path;
    char* output;
    gzFile* handles;
    size_t* nreads;
    FILE* perread;
    FILE* perfile;
} _writer;

typedef _writer* writer;

writer initialize_writer(char* path, char* output_dir, char* perread, char* perfile, char* sample);

void destroy_writer(writer writer);

void write_read(writer writer, kseq_t* seq, size_t barcode);

#endif