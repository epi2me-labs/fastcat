#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include <string.h>

#include "htslib/thread_pool.h"

#include "writer.h"
#include "common.h"
#include "stats.h"
#include "../fastqcomments.h"
#include "bam_utils.h"
#include "../version.h"

// default buffer size for writing to gzip, this is large enough that most reads will not require
// the write buffer to be resized. See gz_fprintf_bridge() below.
// The size is also used with gzbuffer() when opening gzFile handles.
#define GZBUFSIZE 131072  // 128 kB

typedef int (*sink_write_fn)(void *handle, const char *fmt, ...);


static int fprintf_bridge(void *h, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vfprintf((FILE *)h, fmt, ap);
    va_end(ap);
    return r; // <0 on error, otherwise chars printed
}


static int gz_fprintf_bridge(void *h, const char *fmt, ...) {
    size_t bufsize = GZBUFSIZE;
    char *buf = (char*) xalloc(bufsize, sizeof(char), "gzbuffer");

    va_list ap;
    va_start(ap, fmt);
    int need = vsnprintf(buf, bufsize, fmt, ap);
    va_end(ap);

    // only when returned value is non-negative and less than n, the string has been completely written.
    if (need < 0) {
        free(buf);
        return -1;
    }

    if ((size_t)need >= bufsize) {
        bufsize = (size_t) need + 1;
        buf = (char*) xrealloc(buf, bufsize, "gzbuffer");
        va_start(ap, fmt);
        need = vsnprintf(buf, bufsize, fmt, ap);
        va_end(ap);
    }
    int w = gzwrite((gzFile)h, buf, (unsigned)need);
    free(buf);
    return (w == 0) ? -1 : need;  // gzwrite returns <= 0 on error
}


char* strip_path(char* input) {
    if (input == NULL) return NULL;
    size_t len = strlen(input);
    if (input[len - 1] == '/') len--;
    char* output = calloc(len + 1, sizeof(char));
    memcpy(output, input, len);
    return output;
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
     writer->failures = calloc(NUM_FAILURE_CODES, sizeof(uint64_t));
     if (strcmp(sample, "")) {
         // sample is used just for printing to summary, pre-add a tab
         writer->sample = calloc(strlen(sample) + 2, sizeof(char)); 
         strcpy(writer->sample, sample);
         strcat(writer->sample, "\t");
     }
     if (perread != NULL) {
        writer->perread = fopen(perread, "w");
        if (writer->perread == NULL) {
            fprintf(stderr, "Error opening per-read file '%s' for writing.\n", perread);
            exit(EXIT_FAILURE);
        }
        fprintf(writer->perread, "read_id\tfilename\trunid\t");
        if (writer->sample != NULL) fprintf(writer->perread, "sample_name\t");
        fprintf(writer->perread, "read_length\tmean_quality\tchannel\tread_number\tstart_time\n");
     }
     if (perfile != NULL) {
         writer->perfile = fopen(perfile, "w");
         if (writer->perfile == NULL) {
             fprintf(stderr, "Error opening per-file file '%s' for writing.\n", perfile);
             exit(EXIT_FAILURE);
         }
         fprintf(writer->perfile, "filename\t");
         if (writer->sample != NULL) fprintf(writer->perfile, "sample_name\t");
         fprintf(writer->perfile, "n_seqs\tn_bases\tmin_length\tmax_length\tmean_quality");
         for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
             // conver to lowercase for consistency
             const char* src = failure_type[i];
             size_t len = strlen(src);
             char* buffer = calloc(len + 1, sizeof(char));
             for (size_t j = 0; j < len; ++j) {
                 buffer[j] = tolower(src[j]);
             }
             fprintf(writer->perfile, "\t%s", buffer);
             free(buffer);
         }
         fprintf(writer->perfile, "\n");
     }
     if (runids != NULL) {
         writer->runids = fopen(runids, "w");
         if (writer->runids == NULL) {
             fprintf(stderr, "Error opening runids file '%s' for writing.\n", runids);
             exit(EXIT_FAILURE);
         }
         fprintf(writer->runids, "filename\t");
         if (writer->sample != NULL) fprintf(writer->runids, "sample_name\t");
         fprintf(writer->runids, "run_id\tcount\n");
     }
     if (basecallers != NULL) {
         writer->basecallers = fopen(basecallers, "w");
         if (writer->basecallers == NULL) {
             fprintf(stderr, "Error opening basecallers file '%s' for writing.\n", basecallers);
             exit(EXIT_FAILURE);
         }
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
            _write_stats(writer->histograms, writer->output, i, writer->l_stats[i], "length");
            destroy_length_stats(writer->l_stats[i]);
        }

        if(writer->q_stats[i] != NULL) {
            _write_stats(writer->histograms, writer->output, i, writer->q_stats[i], "quality");
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
    free(writer->failures);
    free(writer);
}


void _write_read(writer writer, bam1_t* bam, read_meta meta, void* handle) {
    // choose writing function
    sink_write_fn sink_write = handle == stdout ? fprintf_bridge : gz_fprintf_bridge;

    static const char* fq_comment_fmt = "@%s %s\n%s\n+\n%s\n";
    static const char* sam_comment_fmt = "@%s\t%s\n%s\n+\n%s\n";
    static const char* no_comment_fmt = "@%s\n%s\n+\n%s\n";

    char* qname = bam_get_qname(bam);
    size_t len = bam->core.l_qseq;
    uint8_t* seq_data = bam_get_seq(bam);
    uint8_t* qual_data = bam_get_qual(bam);

    char* seq_str = xalloc(len + 1, sizeof(char), "bam seq");
    char* qual_str = xalloc(len + 1, sizeof(char), "bam qual");
    for (size_t i = 0; i < len; ++i) {
        seq_str[i] = seq_nt16_str[bam_seqi(seq_data, i)];
        qual_str[i] = (qual_data[0] == 0xff) ? '!' : (char)(qual_data[i] + 33);
    }
    seq_str[len] = '\0';
    qual_str[len] = '\0';

    if (writer->reheader) {
        if (meta->tags_str && meta->tags_str->l > 0) {
            sink_write(handle, sam_comment_fmt, qname, meta->tags_str->s, seq_str, qual_str);
        } else {
            sink_write(handle, no_comment_fmt, qname, seq_str, qual_str);
        }
    } else {
        const char* comment = (meta->orig_comment && strlen(meta->orig_comment)) ? meta->orig_comment : meta->comment;
        if (comment && strlen(comment) > 0) {
            sink_write(handle, fq_comment_fmt, qname, comment, seq_str, qual_str);
        } else {
            sink_write(handle, no_comment_fmt, qname, seq_str, qual_str);
        }
    }

    free(seq_str);
    free(qual_str);
}


void _write_read_bam(writer writer, bam1_t* bam, void* handle) {
        if (bam == NULL) {
            fprintf(stderr, "Internal error: BAM record missing for output.\n");
            exit(1);
        }

        if (sam_write1(handle, writer->bam_hdr, bam) < 0) {
            fprintf(stderr, "Error writing read to BAM file.\n");
            exit(1);
        }
}


void _write_stats(char* hist_dir, char* plex_dir, size_t barcode, read_stats* stats, char* type) {
    // write out length stats
    // we assume here the directories have been created already
    char* filepath;
    char* suff = "hist";
    if (plex_dir == NULL) {
        // main output is to stdout, i.e. all read together
        filepath = calloc(strlen(hist_dir) + strlen(type) + strlen(suff) + 3, sizeof(char));
        sprintf(filepath, "%s/%s.%s", hist_dir, type, suff);
    }
    else {
        // demultiplexing
        char* path;
        if (barcode == 0) {
            // unclassified/missing
            path = calloc(strlen(plex_dir) + 15, sizeof(char));
            sprintf(path, "%s/unclassified/", plex_dir);
            filepath = calloc(strlen(path) + strlen(type) + strlen(suff) + 15, sizeof(char));
            sprintf(filepath, "%sunclassified.%s.%s", path, type, suff);
        }
        else {
            path = calloc(strlen(plex_dir) + 14, sizeof(char));
            sprintf(path, "%s/barcode%04lu/", plex_dir, barcode);
            filepath = calloc(strlen(path) + strlen(type) + strlen(suff) + 14, sizeof(char));
            sprintf(filepath, "%sbarcode%04lu.%s.%s", path, barcode, type, suff);
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


void write_read(writer writer, bam1_t* bam, read_meta meta, float mean_q, char* fname) {
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
            bam_get_qname(bam), fname, meta->runid, s, (size_t)bam->core.l_qseq,
                mean_q, meta->channel, meta->read_number, meta->start_time);
    }

    if (writer->output == NULL) {
        // all reads to stdout
        // Stats were initialize in init, no need to check if they exist
        if (writer->write_bam) {
            _write_read_bam(writer, bam, writer->bam_files[0]);
        }
        else {
            _write_read(writer, bam, meta, stdout);
        }
        add_length_count(writer->l_stats[0], bam->core.l_qseq);
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
                _write_read_bam(writer, bam, writer->bam_files[barcode]);
            }
            else {
                // same again for fastq
                if (writer->handles[barcode] == NULL) {
                    create_filepath(writer, barcode, &path, &filepath);
                    ensure_directory(writer, barcode, path);
                    writer->handles[barcode] = gzopen(filepath, "wb");
                    gzbuffer(writer->handles[barcode], GZBUFSIZE);
                }
                _write_read(writer, bam, meta, writer->handles[barcode]);
            }
            free(filepath);
            free(path);
        }

        // handle stats
        if (writer->l_stats[barcode] == NULL) {
            writer->l_stats[barcode] = create_length_stats();
            writer->q_stats[barcode] = create_qual_stats(QUAL_HIST_WIDTH);
        }
        add_length_count(writer->l_stats[barcode], bam->core.l_qseq);
        add_qual_count(writer->q_stats[barcode], mean_q);
        writer->reads_written[barcode]++;
    }
}
