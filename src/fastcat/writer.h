#ifndef FASTCAT_WRITER_H
#define FASTCAT_WRITER_H

#include <zlib.h>


// this gives us kseq_t for below
#ifndef KSEQ_DECLARED
#include "htslib/kseq.h"
KSEQ_DECLARE(gzFile)
#endif

#include "../stats.h"
#include "../fastqcomments.h"

// barcode 0 is reserved for "unclassified"
#define MAX_BARCODES 1025

typedef struct {
    char* path;
    char* output;
    char* histograms;
    gzFile* handles;
    size_t* nreads;
    read_stats** l_stats;
    read_stats** q_stats;
    char* sample;
    size_t reheader;
    FILE* perread;
    FILE* perfile;
    FILE* runids;
    // output file chunking
    size_t reads_per_file;
    size_t* reads_written;
    size_t* file_index;
} _writer;

typedef _writer* writer;

char* strip_path(char* input);

writer initialize_writer(char* output_dir, char* histograms, char* perread, char* perfile, char* runids, char* sample, size_t reheader, size_t reads_per_file);

void destroy_writer(writer writer);

void write_read(writer writer, kseq_t* seq, read_meta meta, float mean_q, char* fname);

#endif
