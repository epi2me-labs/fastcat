#include "regiter.h"
#include "common.h"

int region_from_string(char* input, char** chr, int* start, int* end) {
    *chr = xalloc(strlen(input) + 1, sizeof(char), "chr");
    strcpy(chr, input);
    char *reg_chr = (char *) hts_parse_reg(input, start, end);
    int rtn = 0;
    if (reg_chr) {
        *reg_chr = '\0';  // sets chr to be terminated at correct point
    } else {
        rtn = -1;
    }
    return rtn;
}


int region_from_bed(FILE* bed_fp, char** chr, int* start, int* end) {
    char* line = NULL;
    char* line_copy = NULL;
    size_t len = 0;
    ssize_t read;
    int rtn = 0;

    *start = -1;
    *end = -1;

    if ((read = getline(&line, &len, bed_fp)) != -1) {
        char *newline_pos = strchr(line, '\n');
        if (newline_pos != NULL) {
            *newline_pos = '\0';  // Null-terminate the string at the newline
        }
        line_copy = strdup(line);  // Copy line for error reporting

        // get chromosome
        char* tok = strtok(line, "\t");
        if (tok == NULL) {
            fprintf(stderr, "WARNING: Missing chromosome field in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }
        *chr = xrealloc(*chr, (strlen(tok) + 1) * sizeof(char), "chr");
        strcpy(*chr, tok);

        // get start coordinate
        tok = strtok(NULL, "\t");
        if (tok == NULL) {
            fprintf(stderr, "WARNING: Missing start field in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }
        char* endptr;
        *start = strtol(tok, &endptr, 10);
        if (*endptr != '\0') {
            fprintf(stderr, "WARNING: Invalid start field in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }

        // get end coordinate
        tok = strtok(NULL, "\t");
        if (tok == NULL) {
            fprintf(stderr, "WARNING: Missing end field in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }
        *end = strtol(tok, &endptr, 10);
        if (*endptr != '\0') {
            fprintf(stderr, "WARNING: Invalid end field in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }

        // Validate start and end
        if (*start < 0 || *end < 0 || *end < *start) {
            fprintf(stderr, "WARNING: Invalid region in BED file line: '%s'.\n", line_copy);
            rtn = -2;
            goto cleanup;
        }
    } else {
        rtn = -1;  // EOF
    }

cleanup:
    free(line);
    free(line_copy);
    return rtn;
}



// Initialize the region iterator
regiter init_region_iterator(const char *bed_file, const char *single_region, sam_hdr_t *hdr) {
    regiter it = {0};
    it.hdr = hdr;
    if (bed_file != NULL) {
        it.bed_fp = fopen(bed_file, "r");
        if (it.bed_fp == NULL) {
            fprintf(stderr, "ERROR: Unable to open BED file: %s\n", bed_file);
            it.error |= 1;
        }
    } else if (single_region != NULL) {
        it.single_region = strdup(single_region);
        it.mode = 1; // we'll process this first then switch to BED file mode
    }
    return it;
}


// Clean up the iterator
void destroy_region_iterator(regiter *it) {
    if (it->chr != NULL) free(it->chr);
    if (it->single_region != NULL) free(it->single_region);
    if (it->bed_fp != NULL) fclose(it->bed_fp);
}

// Get the next region
// returns:
//      0 successful
//     -1 if no more regions
//     -2 if error parsing region string
//     -3 if reference not found in BAM header
int next_region(regiter *it) {
    int rtn = 0;
    if (it->mode == 0) {
        rtn = region_from_bed(it->bed_fp, &it->chr, &it->start, &it->end);
    } else if (it->mode == 1) {
        it->mode = 0;
        if (it->single_region) {
            rtn = region_from_string(it->single_region, &it->chr, &it->start, &it->end);
            if (rtn == -1) {
                fprintf(stderr, "WARNING: Failed to parse region string: %s\n", it->single_region);
                it->error |= 2;
                rtn = -2;
            }
        }
    }

    if (rtn == 0) {
        // check reference exists and tidy up length
        int tid = sam_hdr_name2tid(it->hdr, it->chr);
        if (tid < 0) {
            fprintf(stderr, "WARNING: Failed to find reference '%s' in BAM header.\n", it->chr);
            rtn = -3;
        }
        else {
            size_t ref_length = (size_t)sam_hdr_tid2len(it->hdr, tid);
            it->end = min(it->end, (int)ref_length);
            it->n_regions++;
        }
    }
    return rtn;
}
