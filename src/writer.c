#include <sys/stat.h>
#include <sys/types.h>

#include "writer.h"
#include "fastqcomments.h"

char* strip_path(char* input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (input[len - 1] == '/') len--;
    char* output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
}

writer initialize_writer(char* output_dir, char* perread, char* perfile, char* sample, size_t reheader) {
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
     writer->output = strip_path(output_dir);
     writer->handles = calloc(MAX_BARCODES, sizeof(gzFile));
     writer->nreads = calloc(MAX_BARCODES, sizeof(size_t));
     writer->reheader = reheader;
     if (strcmp(sample, "")) {
         // sample is used just for printing to summary, pre-add a tab
         writer->sample = calloc(strlen(sample) + 2, sizeof(char)); 
         strcpy(writer->sample, sample);
         strcat(writer->sample, "\t");
     }
     if (perread != NULL) {
        writer->perread = fopen(perread, "w");
        fprintf(writer->perread, "read_id\tfilename\t");
        if (writer->sample != NULL) fprintf(writer->perread, "sample_name\t");
        fprintf(writer->perread, "read_length\tmean_quality\tchannel\tread_number\tstart_time\n");
     }
     if (perfile != NULL) {
         writer->perfile = fopen(perfile, "w");
         fprintf(writer->perfile, "filename\t");
         if (writer->sample != NULL) fprintf(writer->perfile, "sample_name\t");
         fprintf(writer->perfile, "n_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
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
    if (writer->sample != NULL) free(writer->sample);
    if (writer->perread != NULL) fclose(writer->perread);
    if (writer->perfile != NULL) fclose(writer->perfile);
    if (writer->output != NULL) free(writer->output);
    free(writer->handles);
    free(writer->nreads);
    free(writer);
}

void _write_read(writer writer, kseq_t* seq, read_meta meta, void* handle) {
    int (*write)(void*, const char*, ...) = handle == stdout ? &fprintf : &gzprintf;

    static const char* reheader_fmt = "@%s RD:Z:%s BC:Z:%s BA:Z:%s FC:Z:%s ST:Z:%s RN:i:%lu CH:i:%lu\n%s\n+\n%s\n";
    static const char* wcomment_fmt = "@%s %s\n%s\n+\n%s\n";
    static const char* nocomment_fmt = "@%s\n%s\n+\n%s\n";

    if (writer->reheader && meta->valid) {
        (*write)(
            handle, reheader_fmt,
            seq->name.s, meta->runid, meta->barcode, meta->barcode_alias,
            meta->flow_cell_id, meta->start_time, meta->read_number, meta->channel, 
            seq->seq.s, seq->qual.s);
    }
    else if (seq->comment.l > 0) {
        (*write)(handle, wcomment_fmt, seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
    }
    else {
        (*write)(handle, nocomment_fmt, seq->name.s, seq->seq.s, seq->qual.s);
    }
}


void write_read(writer writer, kseq_t* seq, read_meta meta, float mean_q, char* fname) {
    // TODO: reads per file
    size_t barcode = meta->ibarcode;
    if (barcode > MAX_BARCODES - 1) {
        fprintf(stderr,
            "ERROR: Read's barcode number (%lu) is greater than MAX_BARCODES (%i)\n",
            barcode, MAX_BARCODES);
        exit(1);
    }
    writer->nreads[barcode]++;

    static const char* reheader_fmt = "@%s RD:Z:%s BC:Z:%s BA:Z:%s FC:Z:%s ST:Z:%s RN:i:%lu CH:i:%lu\n%s\n+\n%s\n";
    static const char* wcomment_fmt = "@%s %s\n%s\n+\n%s\n";
    static const char* nocomment_fmt = "@%s\n%s\n+\n%s\n";

    if(writer->perread != NULL) {
        // sample has tab pre-added in init
        char* s = writer->sample == NULL ? "" : writer->sample;
        fprintf(writer->perread, "%s\t%s\t%s%zu\t%1.2f\t%lu\t%lu\t%s\n",
            seq->name.s, fname, s, seq->seq.l, mean_q, meta->channel, meta->read_number, meta->start_time);
    }

    if (writer->output == NULL) {
        _write_read(writer, seq, meta, stdout);
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
        _write_read(writer, seq, meta, writer->handles[barcode]);
    }
}
