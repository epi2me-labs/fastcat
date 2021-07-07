#include <sys/stat.h>
#include <sys/types.h>

#include "writer.h"

writer initialize_writer(char* path, char* output_dir, char* perread, char* perfile, char* sample) {
    if (output_dir != NULL) {
        int rtn = mkdir(output_dir, 0700);
        if (rtn == -1) {
            fprintf(stderr,
               "Error: Cannot create output directory '%s'. Check location is writeable and does not already exist.\n",
               output_dir);
            return NULL;
        }
     }
     writer writer = calloc(1, sizeof(_writer));
     writer->output = output_dir; // TODO: take a copy
     writer->handles = calloc(MAX_BARCODES, sizeof(gzFile));
     writer->nreads = calloc(MAX_BARCODES, sizeof(size_t));
     if (perread != NULL) {
         writer->perread = fopen(perread, "w");
        if (strcmp(sample, "")) {
            fprintf(writer->perread, "read_id\tfilename\tsample_name\tread_length\tmean_quality\n");
        } else {
            fprintf(writer->perread, "read_id\tfilename\tread_length\tmean_quality\n");
        }
     }
     if (perfile != NULL) {
         writer->perfile = fopen(perfile, "w");
         if (strcmp(sample, "")) {
             fprintf(writer->perfile, "filename\tsample_name\tn_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
         } else {
             fprintf(writer->perfile, "filename\tn_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
         }
     }
     return writer;
}


void destroy_writer(writer writer) {
    for(size_t i=0; i < MAX_BARCODES; ++i) {
       if(writer->handles[i] != NULL) {
           gzflush(writer->handles[i], Z_FINISH);
           gzclose(writer->handles[i]);
       }
    }
    if (writer->perread != NULL) {
        fclose(writer->perread);
    }
    if (writer->perfile != NULL) {
        fclose(writer->perfile);
    }
    free(writer->handles);
    free(writer->nreads);
    free(writer);
}


void write_read(writer writer, kseq_t* seq, size_t barcode) {
    // TODO: reads per file
    if (barcode > MAX_BARCODES - 1) {
        fprintf(stderr,
            "ERROR: Read's barcode number (%lu) is greater than MAX_BARCODES (%i)\n",
            barcode, MAX_BARCODES);
        exit(1);
    }
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
