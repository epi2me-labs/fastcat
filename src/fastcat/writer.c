#include <sys/stat.h>
#include <sys/types.h>

#include "writer.h"
#include "common.h"
#include "stats.h"
#include "../fastqcomments.h"

// default buffer size for writing to gzip, this is large enough that most reads will not require
// the write buffer to be resized. See _gzsnprintf() below.
// The size is also used with gzbuffer() when opening gzFile handles.
#define GZBUFSIZE 131072  // 128 kB



char* strip_path(char* input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (input[len - 1] == '/') len--;
    char* output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
}


// Safely write a formatted string to a gzFile pointer
int _gzsnprintf(gzFile file, const char *format, ...) {
    va_list myargs;
    va_start(myargs, format);
    int bufsize = GZBUFSIZE;
    char* buf = (char*) xalloc(bufsize, sizeof(char), "gzbuffer");
    int written = vsnprintf(buf, bufsize, format, myargs);
    va_end(myargs);
    if (written > bufsize) {
        bufsize = written;
        buf = (char *) realloc(buf, bufsize);
        va_start(myargs, format);
        written = vsnprintf(buf, bufsize, format, myargs);
        va_end(myargs);
    }
    gzputs(file, buf);
    free(buf);
    return written;
}


writer initialize_writer(char* output_dir, char* histograms, char* perread, char* perfile, char* runids, char* sample, size_t reheader, size_t reads_per_file) {
    if (output_dir != NULL) {  // demultiplexing
        int rtn = mkdir_hier(output_dir);
        if (rtn == -1) {
            fprintf(stderr,
               "Error: Cannot create output directory '%s'. Check location is writeable and directory does not exist.\n",
               output_dir);
            return NULL;
        }
     }
     else {
         // histograms go in their own directory when not demultiplexing
         int rtn = mkdir_hier(histograms);
         if (rtn == -1) {
            fprintf(stderr,
                "Error: Cannot create output directory '%s'. Check location is writeable and directory does not exist.\n",
                histograms);
            exit(EXIT_FAILURE);
         }
     }

     writer writer = calloc(1, sizeof(_writer));
     writer->output = strip_path(output_dir);
     writer->histograms = strip_path(histograms);
     writer->handles = calloc(MAX_BARCODES, sizeof(gzFile));
     writer->nreads = calloc(MAX_BARCODES, sizeof(size_t));
     writer->l_stats = calloc(MAX_BARCODES, sizeof(read_stats*));
     writer->q_stats = calloc(MAX_BARCODES, sizeof(read_stats*));
     // we want at least 1 stats accumulator (when not demultiplexing)
     writer->l_stats[0] = create_length_stats();
     writer->q_stats[0] = create_qual_stats(QUAL_HIST_WIDTH);
     writer->reheader = reheader;
     writer->reads_per_file = reads_per_file;
     writer->reads_written = calloc(MAX_BARCODES, sizeof(size_t));
     writer->file_index = calloc(MAX_BARCODES, sizeof(size_t));
     if (strcmp(sample, "")) {
         // sample is used just for printing to summary, pre-add a tab
         writer->sample = calloc(strlen(sample) + 2, sizeof(char)); 
         strcpy(writer->sample, sample);
         strcat(writer->sample, "\t");
     }
     if (perread != NULL) {
        writer->perread = fopen(perread, "w");
        fprintf(writer->perread, "read_id\tfilename\trunid\t");
        if (writer->sample != NULL) fprintf(writer->perread, "sample_name\t");
        fprintf(writer->perread, "read_length\tmean_quality\tchannel\tread_number\tstart_time\n");
     }
     if (perfile != NULL) {
         writer->perfile = fopen(perfile, "w");
         fprintf(writer->perfile, "filename\t");
         if (writer->sample != NULL) fprintf(writer->perfile, "sample_name\t");
         fprintf(writer->perfile, "n_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
     }
     if (runids != NULL) {
         writer->runids = fopen(runids, "w");
         fprintf(writer->runids, "filename\t");
         if (writer->sample != NULL) fprintf(writer->runids, "sample_name\t");
         fprintf(writer->runids, "run_id\tcount\n");
     }
     return writer;
}


void _write_stats(char* hist_dir, char* plex_dir, size_t barcode, read_stats* stats, char* type);

void destroy_writer(writer writer) {
    for(size_t i=0; i < MAX_BARCODES; ++i) {
        if(writer->handles[i] != NULL) {
            gzflush(writer->handles[i], Z_FINISH);
            gzclose(writer->handles[i]);
        }

        if(writer->l_stats[i] != NULL) {
            _write_stats(writer->histograms, writer->output, i, writer->l_stats[i], "length\0");
            destroy_length_stats(writer->l_stats[i]);
        }

        if(writer->q_stats[i] != NULL) {
            _write_stats(writer->histograms, writer->output, i, writer->q_stats[i], "quality\0");
            destroy_qual_stats(writer->q_stats[i]);
        }

    }
    if (writer->sample != NULL) free(writer->sample);
    if (writer->perread != NULL) fclose(writer->perread);
    if (writer->perfile != NULL) fclose(writer->perfile);
    if (writer->runids != NULL) fclose(writer->runids);
    if (writer->output != NULL) free(writer->output);
    if (writer->histograms != NULL) free(writer->histograms);
    free(writer->handles);
    free(writer->nreads);
    free(writer->l_stats);
    free(writer->q_stats);
    free(writer->reads_written);
    free(writer->file_index);
    free(writer);
}

void _write_read(writer writer, kseq_t* seq, read_meta meta, void* handle) {
    
    int (*write)(void*, const char*, ...);
    if (handle == stdout) { write = &fprintf; } else { write = &_gzsnprintf; }

    static const char* fq_comment_fmt = "@%s %s\n%s\n+\n%s\n";
    static const char* sam_comment_fmt = "@%s\t%s\n%s\n+\n%s\n";
    static const char* no_comment_fmt = "@%s\n%s\n+\n%s\n";

    if (seq->comment.l > 0) {
        if (writer->reheader) {
            (*write)(handle, sam_comment_fmt, seq->name.s, meta->tags_str->s, seq->seq.s, seq->qual.s);
        }
        else {
            (*write)(handle, fq_comment_fmt, seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
        }
    }
    else {
        (*write)(handle, no_comment_fmt, seq->name.s, seq->seq.s, seq->qual.s);
    }
}


void _write_stats(char* hist_dir, char* plex_dir, size_t barcode, read_stats* stats, char* type) {
    // write out length stats
    // we assume here the directories have been created already
    char* filepath;
    if (plex_dir == NULL) {
        // main output is to stdout, i.e. all read together
        filepath = calloc(strlen(hist_dir) + strlen(type) + 7, sizeof(char));
        sprintf(filepath, "%s/%s.hist", hist_dir, type);
    }
    else {
        // demultiplexing
        char* path;
        if (barcode == 0) {
            // unclassified/missing
            path = calloc(strlen(plex_dir) + 15, sizeof(char));
            sprintf(path, "%s/unclassified/", plex_dir);
            filepath = calloc(strlen(path) + strlen(type) + 19, sizeof(char));
            sprintf(filepath, "%sunclassified.%s.hist", path, type);
        }
        else {
            path = calloc(strlen(plex_dir) + 14, sizeof(char));
            sprintf(path, "%s/barcode%04lu/", plex_dir, barcode);
            filepath = calloc(strlen(path) + strlen(type) + 18, sizeof(char));
            sprintf(filepath, "%sbarcode%04lu.%s.hist", path, barcode, type);
        }
        free(path);
    }
    FILE* fp = fopen(filepath, "w");
    print_stats(stats, false, true, fp);
    fclose(fp);
    free(filepath);
}


void create_filepath(writer writer, size_t barcode, char** path, char** filepath) {
    // additional file index string if needed
    // we'll allocate just an empty string
    int ex_size = (writer->reads_per_file == 0) ? 0 : 5;
    char* extra = calloc(ex_size + 1, sizeof(char));
    if (writer->reads_per_file != 0) {
        sprintf(extra, "_%04zu", writer->file_index[barcode]);
    }

    if (barcode == 0) { // unclassified/missing
        *path = (char*)calloc(strlen(writer->output) + 15, sizeof(char));
        sprintf(*path, "%s/unclassified/", writer->output);
        *filepath = (char*)calloc(strlen(*path) + 22 + ex_size, sizeof(char));
        sprintf(*filepath, "%sunclassified%s.fastq.gz", *path, extra);
    }
    else { // barcoded data
        *path = (char*)calloc(strlen(writer->output) + 14, sizeof(char));
        sprintf(*path, "%s/barcode%04lu/", writer->output, barcode);
        *filepath = (char*)calloc(strlen(*path) + 21 + ex_size, sizeof(char));
        sprintf(*filepath, "%sbarcode%04lu%s.fastq.gz", *path, barcode, extra);
    }
    free(extra);
}


void write_read(writer writer, kseq_t* seq, read_meta meta, float mean_q, char* fname) {
    size_t barcode = meta->ibarcode;
    if (barcode > MAX_BARCODES - 1) {
        fprintf(stderr,
            "ERROR: Read's barcode number (%lu) is greater than MAX_BARCODES (%i)\n",
            barcode, MAX_BARCODES);
        exit(1);
    }
    writer->nreads[barcode]++;

    if(writer->perread != NULL) {
        // sample has tab pre-added in init
        char* s = writer->sample == NULL ? "" : writer->sample;
        fprintf(writer->perread, "%s\t%s\t%s\t%s%zu\t%.2f\t%lu\t%lu\t%s\n",
            seq->name.s, fname, meta->runid, s, seq->seq.l, \
                mean_q, meta->channel, meta->read_number, meta->start_time);
    }

    if (writer->output == NULL) {
        // all reads to stdout
        _write_read(writer, seq, meta, stdout);
        add_length_count(writer->l_stats[0], seq->seq.l);
        add_qual_count(writer->q_stats[0], mean_q);
    }
    else {
        // demultiplexing reads
        // first handle multipart-output 
        if (writer->reads_per_file != 0 && writer->reads_written[barcode] == writer->reads_per_file) {
            if (writer->handles[barcode] == NULL) {
                fprintf(stderr, "Unexpected output file status encountered.");
                exit(1);
            }
            gzflush(writer->handles[barcode], Z_FINISH);
            gzclose(writer->handles[barcode]);
            writer->handles[barcode] = NULL;
            writer->file_index[barcode]++;
            writer->reads_written[barcode] = 0;
        }

        // open a file, if we need to
        if (writer->handles[barcode] == NULL) {
            char* path = NULL;
            char* filepath = NULL;
            create_filepath(writer, barcode, &path, &filepath);
           
            // if single file, or first file and no reads yet, make the directory
            if (writer->reads_per_file == 0 
                    || (writer->file_index[barcode] == 0 && writer->reads_written[barcode] == 0)) {
                int rtn = mkdir_hier(path);
                if (rtn == -1) {
                    fprintf(stderr, "Failed to create barcode directory '%s\n'.", path);
                    exit(1);
                }
            }
            writer->handles[barcode] = gzopen(filepath, "wb");
            gzbuffer(writer->handles[barcode], GZBUFSIZE);
            free(filepath);
            free(path);
        }

        // all of the above was just to get a handle (on life), now use it
        _write_read(writer, seq, meta, writer->handles[barcode]);

        // handle stats
        if (writer->l_stats[barcode] == NULL) {
            writer->l_stats[barcode] = create_length_stats();
            writer->q_stats[barcode] = create_qual_stats(QUAL_HIST_WIDTH);
        }
        add_length_count(writer->l_stats[barcode], seq->seq.l);
        add_qual_count(writer->q_stats[barcode], mean_q);
        writer->reads_written[barcode]++;
    }
}
