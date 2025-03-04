#ifndef FASTCAT_WRITER_H
#define FASTCAT_WRITER_H

#include <zlib.h>


// this gives us kseq_t for below
#ifndef KSEQ_DECLARED
#include "htslib/kseq.h"
KSEQ_DECLARE(gzFile)
#endif

#include <htslib/sam.h> // HTSlib for BAM output

#include "../stats.h"
#include "../fastqcomments.h"

// barcode 0 is reserved for "unclassified"
#define MAX_BARCODES 1025

typedef struct {
    char* output;
    char* histograms;
    gzFile* handles;
    size_t* nreads;
    size_t* reads_written;
    size_t* file_index;
    read_stats** l_stats;
    read_stats** q_stats;
    FILE* perread;
    FILE* perfile;
    FILE* runids;
    FILE* basecallers;
    char* sample;
    size_t reheader;
    size_t reads_per_file;
    // optional BAM conversion
    int write_bam;
    htsFile** bam_files;
    bam_hdr_t* bam_hdr;
    htsThreadPool hts_pool;
} _writer;

typedef _writer* writer;

char* strip_path(char* input);

writer initialize_writer(
        char* output_dir, char* histograms, char* perread, char* perfile,
        char* runids, char* basecallers, char* sample,
        size_t reheader, size_t write_bam, size_t reads_per_file,
        int threads);

void destroy_writer(writer writer);

void write_read(writer writer, kseq_t* seq, read_meta meta, float mean_q, char* fname);

#endif
