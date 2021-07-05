
#include "writer.h"

writer initialize_writer(char* path, size_t to_stdout) {
     writer writer = calloc(1, sizeof(_writer));
     writer->to_stdout = to_stdout;
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


void write_read(writer writer, kseq_t* seq, size_t barcode, char* alias) {
    // TODO: reads per file
    writer->nreads[barcode]++;
    if (writer->to_stdout) {
        if (seq->comment.l > 0) {
            fprintf(stdout, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
        }
        else {
            fprintf(stdout, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
        }
    }
    else {
        if (writer->handles[barcode] == NULL) {
            char* filepath = calloc(21, sizeof(char));
            sprintf(filepath, "barcode%04lu.fastq.gz", barcode);
            writer->handles[barcode] = gzopen(filepath, "wb");
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
