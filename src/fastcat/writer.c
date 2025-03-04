#include <sys/stat.h>
#include <sys/types.h>

#include "htslib/thread_pool.h"

#include "writer.h"
#include "common.h"
#include "stats.h"
#include "../fastqcomments.h"
#include "../version.h"

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


writer initialize_writer(
        char* output_dir, char* histograms, char* perread, char* perfile,
        char* runids, char* basecallers, char* sample,
        size_t reheader, size_t write_bam, size_t reads_per_file,
        int threads) {
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

     writer->nreads = calloc(MAX_BARCODES, sizeof(size_t));
     writer->l_stats = calloc(MAX_BARCODES, sizeof(read_stats*));
     writer->q_stats = calloc(MAX_BARCODES, sizeof(read_stats*));
     // we want a stats accumulator (when not demultiplexing)
     // This also ensures we write out an empty histogram file when
     // no reads are processed. (To go with our other empty summary files)
     if (writer->output == NULL) {
         writer->l_stats[0] = create_length_stats();
         writer->q_stats[0] = create_qual_stats(QUAL_HIST_WIDTH);
     }
     writer->reheader = reheader;
     writer->write_bam = write_bam;
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
     if (basecallers != NULL) {
         writer->basecallers = fopen(basecallers, "w");
         fprintf(writer->basecallers, "filename\t");
         if (writer->sample != NULL) fprintf(writer->basecallers, "sample_name\t");
         fprintf(writer->basecallers, "basecaller\tcount\n");
     }

     // we alloc all the file pointers here, but we might not use them, just to keep the code simple
     // whats MAX_BARCODES * a few bytes between friends?
     if (write_bam) {
             fprintf(stderr, "Using %d threads for BAM writing\n", threads);
             writer->hts_pool.pool = hts_tpool_init(threads);
             writer->hts_pool.qsize = 0;
             if (writer->hts_pool.pool == NULL) {
                 fprintf(stderr, "Error creating thread pool\n");
                 exit(1);
             }
             // later...call hts_set_opt on each fp opened

         writer->bam_hdr = sam_hdr_init();
         sam_hdr_add_line(writer->bam_hdr, "HD", "VN", SAM_FORMAT_VERSION, "SO", "unsorted", NULL);
         sam_hdr_add_line(writer->bam_hdr, "PG", "ID", "fastcat", "PN", "fastcat", "VN", argp_program_version, NULL);
         writer->bam_files = calloc(MAX_BARCODES, sizeof(htsFile*));
         if (writer->output == NULL) { // to stdout
             writer->bam_files[0] = hts_open("-", "wb");
             hts_set_opt(writer->bam_files[0], HTS_OPT_THREAD_POOL, &writer->hts_pool);
             if (sam_hdr_write(writer->bam_files[0], writer->bam_hdr)) {
                 fprintf(stderr, "Error writing header to BAM on stdout\n");
                 exit(1);
             }
         }
     }
     else { // fastq output
        writer->handles = calloc(MAX_BARCODES, sizeof(gzFile));
     }

     return writer;
}


void _write_stats(char* hist_dir, char* plex_dir, size_t barcode, read_stats* stats, char* type);


void destroy_writer(writer writer) {
    for(size_t i=0; i < MAX_BARCODES; ++i) {
        if (writer->write_bam) {
            if (writer->bam_files[i] != NULL) {
                hts_close(writer->bam_files[i]);
            }
        }
        else {
            if (writer->handles[i] != NULL) {
                gzflush(writer->handles[i], Z_FINISH);
                gzclose(writer->handles[i]);
            }
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
    if (writer->write_bam) { // must be after file closing
        hts_tpool_destroy(writer->hts_pool.pool);
    }

    if (writer->sample != NULL) free(writer->sample);
    if (writer->perread != NULL) fclose(writer->perread);
    if (writer->perfile != NULL) fclose(writer->perfile);
    if (writer->runids != NULL) fclose(writer->runids);
    if (writer->basecallers != NULL) fclose(writer->basecallers);
    if (writer->output != NULL) free(writer->output);
    if (writer->histograms != NULL) free(writer->histograms);
    if (writer->write_bam) {
        free(writer->bam_files);
        bam_hdr_destroy(writer->bam_hdr);
    }
    else {
        free(writer->handles);
    }
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


// htslib has aux_parse but its static :(
int parse_and_set_aux_tags(bam1_t *b, const char *aux_str) {
    if (!b || !aux_str) {
        fprintf(stderr, "Invalid input to parse_and_set_aux_tags\n");
        return -1;
    }

    const char *ptr = aux_str;
    while (*ptr) {
        char tag[3] = {0};
        char type;
        int consumed = 0;

        // Read the tag and type
        if (sscanf(ptr, "%2s:%c%n", tag, &type, &consumed) != 2) {
            fprintf(stderr, "Invalid aux format: %s\n", ptr);
            return -1;
        }
        ptr += consumed+1;

        size_t len = 0;
        int int_value;
        float float_value;
        if (type == 'i' || type == 'I' || type == 'c' || type == 'C' || type == 's' || type == 'S') {
            char *endptr;
            int_value = strtol(ptr, &endptr, 10);
            if (endptr == ptr) {
                fprintf(stderr, "Failed to parse integer for tag %s\n", tag);
                return -1;
            }
            bam_aux_update_int(b, tag, int_value);
            ptr = endptr;
        }
        else if (type == 'f') {
            char *endptr;
            float_value = strtof(ptr, &endptr);
            if (endptr == ptr) {
                fprintf(stderr, "Failed to parse float for tag %s\n", tag);
                return -1;
            }
            bam_aux_update_float(b, tag, float_value);
            ptr = endptr;
        }
        else if (type == 'Z' || type == 'H') {
            len = strcspn(ptr, "\t\n");
            bam_aux_update_str(b, tag, len, ptr);
            ptr += len;
        }
        else {
            fprintf(stderr, "Unsupported tag type: %c\n", type);
            return -1;
        }

        while (*ptr == '\t' || *ptr == ' ') {
            ptr++;
        }
    }
    
    return 0;
}


void _write_read_bam(writer writer, kseq_t* seq, read_meta meta, void* handle) {
        bam1_t* b = bam_init1();
        bam_set1(
            b,
            seq->name.l, seq->name.s,
            4, -1, -1, 0,
            0, NULL,
            -1, -1, 0,
            seq->seq.l, seq->seq.s, seq->qual.s,
            0
        );
        // the implementation of bam_set1() seems not to take into account the 33 offset
        // you'd typically have in a string encoding
        uint8_t* qual = b->data + b->core.l_qname + (b->core.l_qseq + 1) / 2;
        for (int i = 0; i < b->core.l_qseq; ++i) {
            qual[i] -= 33;
        }

        // see fastqcomments.c for the definition of read_meta
        // there we parcelled everything up nicely into a sam formatted string
        // with garbage being dumped into CO:Z tag
        // setting the known tags directly is futile as they would be overwritten
        // by this function in any case
        if (meta->tags_str->l > 0) {
            if (parse_and_set_aux_tags(b, meta->tags_str->s) < 0) {
                fprintf(stderr, "Error parsing auxiliary tags\n");
                fprintf(stderr, "read: %s\n", seq->name.s);
                fprintf(stderr, "tags: %s\n", meta->tags_str->s);
                fprintf(stderr, "rest: %s\n", meta->rest->s);
                exit(1);
            }
        }

        if (sam_write1(handle, writer->bam_hdr, b) < 0) {
            fprintf(stderr, "Error writing read to BAM file.\n");
            exit(1);
        }
        bam_destroy1(b);
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
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filepath);
        exit(1);
    }
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
        sprintf(*filepath, "%sunclassified%s.%s", *path, extra, writer->write_bam ? "bam" : "fastq.gz");
    }
    else { // barcoded data
        *path = (char*)calloc(strlen(writer->output) + 14, sizeof(char));
        sprintf(*path, "%s/barcode%04lu/", writer->output, barcode);
        *filepath = (char*)calloc(strlen(*path) + 21 + ex_size, sizeof(char));
        sprintf(*filepath, "%sbarcode%04lu%s.%s", *path, barcode, extra, writer->write_bam ? "bam" : "fastq.gz");
    }
    free(extra);
}


void ensure_directory(writer writer, size_t barcode, char* path) {
    // if single file, or first file and no reads yet, make the directory
    if (writer->reads_per_file == 0 
            || (writer->file_index[barcode] == 0 && writer->reads_written[barcode] == 0)) {
        int rtn = mkdir_hier(path);
        if (rtn == -1) {
            fprintf(stderr, "Failed to create barcode directory '%s\n'.", path);
            exit(1);
        }
    }
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
        // Stats were initialize in init, no need to check if they exist
        if (writer->write_bam) {
            _write_read_bam(writer, seq, meta, writer->bam_files[0]);
        }
        else {
            _write_read(writer, seq, meta, stdout);
        }
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
            if (writer->write_bam) {
                hts_close(writer->bam_files[barcode]);
                writer->bam_files[barcode] = NULL;
            }
            else {
                gzflush(writer->handles[barcode], Z_FINISH);
                gzclose(writer->handles[barcode]);
                writer->handles[barcode] = NULL;
            }
            writer->file_index[barcode]++;
            writer->reads_written[barcode] = 0;
        }

        // write read to correct file
        {
            char* path = NULL;
            char* filepath = NULL;
            if (writer->write_bam) {
                // open a file, if we need to
                if (writer->bam_files[barcode] == NULL) {
                    create_filepath(writer, barcode, &path, &filepath);
                    ensure_directory(writer, barcode, path);
                    writer->bam_files[barcode] = hts_open(filepath, "wb");
                    hts_set_opt(writer->bam_files[barcode], HTS_OPT_THREAD_POOL, &writer->hts_pool);
                    if (sam_hdr_write(writer->bam_files[barcode], writer->bam_hdr)) {
                        fprintf(stderr, "Error writing header to BAM file\n");
                        exit(1);
                    }
                }
                _write_read_bam(writer, seq, meta, writer->bam_files[barcode]);
            }
            else {
                // same again for fastq
                if (writer->handles[barcode] == NULL) {
                    create_filepath(writer, barcode, &path, &filepath);
                    ensure_directory(writer, barcode, path);
                    writer->handles[barcode] = gzopen(filepath, "wb");
                    gzbuffer(writer->handles[barcode], GZBUFSIZE);
                }
                _write_read(writer, seq, meta, writer->handles[barcode]);
            }
            free(filepath);
            free(path);
        }

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
