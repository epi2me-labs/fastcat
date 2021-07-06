#include <sys/stat.h>
#include <sys/types.h>

#include "writer.h"

writer initialize_writer(char* path, char* output_dir) {
     writer writer = calloc(1, sizeof(_writer));
     writer->output = output_dir; // TODO: take a copy
     writer->handles = calloc(MAX_BARCODES, sizeof(gzFile));
     writer->nreads = calloc(MAX_BARCODES, sizeof(size_t));
     return writer;
}


void destroy_writer(writer writer) {
    for(size_t i=0; i < MAX_BARCODES; ++i) {
       if(writer->handles[i] != NULL) {
           gzflush(writer->handles[i], Z_FINISH);
           gzclose(writer->handles[i]);
       }
    }
    free(writer->handles);
    free(writer->nreads);
    free(writer);
}


void write_read(writer writer, kseq_t* seq, size_t barcode) {
    // TODO: reads per file
    writer->nreads[barcode]++;
    if (writer->output == NULL) {
        if (seq->comment.l > 0) {
            fprintf(stdout, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
        }
        else {
            fprintf(stdout, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
        }
    }
    else {
        if (writer->handles[barcode] == NULL) {
            char* path;
            char* filepath;
            if (barcode == 0) {
                // unclassified/missing
                path = calloc(strlen(writer->output) + 15, sizeof(char));
                sprintf(path, "%s/unclassified/", writer->output);
                filepath = calloc(strlen(path) + 22, sizeof(char));
                sprintf(filepath, "%sunclassified.fastq.gz", path);
            }
            else {
                path = calloc(strlen(writer->output) + 14, sizeof(char));
                sprintf(path, "%s/barcode%04lu/", writer->output, barcode);
                filepath = calloc(strlen(path) + 21, sizeof(char));
                sprintf(filepath, "%sbarcode%04lu.fastq.gz", path, barcode);
            }
            
            int rtn = mkdir(path, 0700);
            if (rtn == -1) {
                fprintf(stderr, "Failed to create barcode directory '%s\n'.", path);
                exit(1);
            }
            writer->handles[barcode] = gzopen(filepath, "wb");
            free(path);
            free(filepath);
        }
        gzFile handle = writer->handles[barcode];
        ssize_t nbytes;
        if (seq->comment.l > 0) {
            gzprintf(writer->handles[barcode], "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
        }
        else {
            gzprintf(writer->handles[barcode], "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
        }
    }
}
