#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kseq.h"
#include <argp.h>

KSEQ_INIT(gzFile, gzread)

static inline size_t max ( size_t a, size_t b ) { return a > b ? a : b; }
static inline size_t min ( size_t a, size_t b ) { return a < b ? a : b; }

const float qprobs[100] = {
	1.00000000e+00, 7.94328235e-01, 6.30957344e-01, 5.01187234e-01,
    3.98107171e-01, 3.16227766e-01, 2.51188643e-01, 1.99526231e-01,
    1.58489319e-01, 1.25892541e-01, 1.00000000e-01, 7.94328235e-02,
    6.30957344e-02, 5.01187234e-02, 3.98107171e-02, 3.16227766e-02,
    2.51188643e-02, 1.99526231e-02, 1.58489319e-02, 1.25892541e-02,
    1.00000000e-02, 7.94328235e-03, 6.30957344e-03, 5.01187234e-03,
    3.98107171e-03, 3.16227766e-03, 2.51188643e-03, 1.99526231e-03,
    1.58489319e-03, 1.25892541e-03, 1.00000000e-03, 7.94328235e-04,
    6.30957344e-04, 5.01187234e-04, 3.98107171e-04, 3.16227766e-04,
    2.51188643e-04, 1.99526231e-04, 1.58489319e-04, 1.25892541e-04,
    1.00000000e-04, 7.94328235e-05, 6.30957344e-05, 5.01187234e-05,
    3.98107171e-05, 3.16227766e-05, 2.51188643e-05, 1.99526231e-05,
    1.58489319e-05, 1.25892541e-05, 1.00000000e-05, 7.94328235e-06,
    6.30957344e-06, 5.01187234e-06, 3.98107171e-06, 3.16227766e-06,
    2.51188643e-06, 1.99526231e-06, 1.58489319e-06, 1.25892541e-06,
    1.00000000e-06, 7.94328235e-07, 6.30957344e-07, 5.01187234e-07,
    3.98107171e-07, 3.16227766e-07, 2.51188643e-07, 1.99526231e-07,
    1.58489319e-07, 1.25892541e-07, 1.00000000e-07, 7.94328235e-08,
    6.30957344e-08, 5.01187234e-08, 3.98107171e-08, 3.16227766e-08,
    2.51188643e-08, 1.99526231e-08, 1.58489319e-08, 1.25892541e-08,
    1.00000000e-08, 7.94328235e-09, 6.30957344e-09, 5.01187234e-09,
    3.98107171e-09, 3.16227766e-09, 2.51188643e-09, 1.99526231e-09,
    1.58489319e-09, 1.25892541e-09, 1.00000000e-09, 7.94328235e-10,
    6.30957344e-10, 5.01187234e-10, 3.98107171e-10, 3.16227766e-10,
    2.51188643e-10, 1.99526231e-10, 1.58489319e-10, 1.25892541e-10};


float mean_qual(char* qual, size_t len) {
	float qsum = 0;
	for (size_t i=0; i<len; ++i) {
		//printf("%s\n", qual);
		int q = (int)(qual[i]) - 33;
		//printf("%c\t%i\n", q, qual);
		qsum += qprobs[(int) (qual[i]) - 33];
	}
    qsum /= len;
	return -10 * log10(qsum);
}

const char *argp_program_version = "0.1.0";
const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
  "fastcat -- concatenate and summarise .fastq(.gz) files.\
  \vInput files may be given on stdin by specifing the input as '-'.";
static char args_doc[] = "reads1.fastq(.gz) reads2.fastq(.gz) ...";
static struct argp_option options[] = {
    {"read",   'r',  "READ SUMMARY",  0,  "Per-read summary output"},
    {"file",   'f',  "FILE SUMMARY",  0,  "Per-file summary output"},
    {"sample", 's',  "SAMPLE NAME",   0,  "Sample name (if given adds a 'sample_name' column)"},
    { 0 }
};

struct arguments {
    char *perread;
    char *perfile;
    char *sample;
    char **files;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch (key) {
        case 's':
            arguments->sample = arg;
            break;
        case 'r':
            arguments->perread = arg;
            break;
        case 'f':
            arguments->perfile = arg;
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            arguments->files = &state->argv[state->next - 1];
            state->next = state->argc;
            break;
        defualt:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};


int process_file(char* fname, char* sample, FILE* outfp, FILE* summaryfp) {

	gzFile fp;
	kseq_t *seq;
        
    fp = gzopen(fname, "r");
    seq = kseq_init(fp);
    size_t n = 0, slen=0;
    size_t minl=UINTMAX_MAX, maxl=0;
    double meanq = 0; // this may lose precision
    while (kseq_read(seq) >= 0) {
        ++n ; slen += seq->seq.l;
        minl = min(minl, seq->seq.l);
        maxl = max(maxl, seq->seq.l);
    	if (seq->qual.l == 0) {
    		fprintf(stderr, "No quality string found for '%s' (FASTA is unsupported).\n", seq->name.s);
    		return 1;
    	}
    	float mean_q = mean_qual(seq->qual.s, seq->qual.l);
        meanq += mean_q;
        fprintf(outfp, "%s\t%s\t%s%zu\t%1.2f\n", seq->name.s, fname, sample, seq->seq.l, mean_q);
    	if (seq->comment.l > 0) {
    		fprintf(stdout, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
    	} else {
    		fprintf(stdout, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
    	}
    }
    fprintf(summaryfp, "%s\t%s%zu\t%zu\t%zu\t%zu\t%1.2f\n", fname, sample, n, slen, minl, maxl, meanq/n);
    kseq_destroy(seq);
    gzclose(fp);
}


int main(int argc, char **argv) {
    struct arguments args;
    args.perread = "read-summary.txt";
    args.perfile = "file-summary.txt";
    args.sample = "";

    argp_parse(&argp, argc, argv, 0, 0, &args);
    char *sample;
    if (strcmp(args.sample, "")) {
        fprintf(stderr, "Adding sample\n");
        sample = calloc(strlen(args.sample) + 2, sizeof(char)); 
        strcpy(sample, args.sample);
        strcat(sample, "\t");
    } else {
        sample = "";
    }

    int nfile = 0;
    for( ; args.files[nfile] ; nfile++);

    FILE* outfp = fopen(args.perread, "w");
    FILE* summaryfp = fopen(args.perfile, "w");
    if (strcmp(args.sample, "")) {
        fprintf(outfp, "read_id\tfilename\tsample_name\tread_length\tmean_quality\n");
	    fprintf(summaryfp, "filename\tsample_name\tn_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
    } else {
        fprintf(outfp, "read_id\tfilename\tread_length\tmean_quality\n");
	    fprintf(summaryfp, "filename\tn_seqs\tn_bases\tmin_length\tmax_length\tmean_quality\n");
    }

    if (nfile==1 && strcmp(args.files[0], "-") == 0) {
        char *ln = NULL;
        size_t n = 0;
        ssize_t nchr = 0;
        while ((nchr = getline (&ln, &n, stdin)) != -1) {
            ln[strcspn(ln, "\r\n")] = 0;
            int rtn = process_file(ln, sample, outfp, summaryfp);
            if (rtn != 0) return rtn;
        }
        free(ln);
    } else { 
        for (size_t i=0; i<nfile; ++i) {
            int rtn = process_file(args.files[i], sample, outfp, summaryfp);
            if (rtn != 0) return rtn;
        }
    }
	fclose(outfp);
	fclose(summaryfp);

	return 0;
}
