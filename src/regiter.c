#include "regiter.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "fc_log.h"

int region_from_string(char *input, char **chr, int *start, int *end)
{
    int rtn = 0;
    char *tmp = fc_strdup(input, "region string for parsing");
    char *reg_chr = (char *)hts_parse_reg(tmp, start, end);
    if (reg_chr) {
        *reg_chr = '\0'; // terminate at region separator
        size_t len = strlen(tmp);
        *chr = fc_alloc(len + 1, sizeof(char), "chr");
        memcpy(*chr, tmp, len + 1);
    } else {
        rtn = -1;
    }
    free(tmp);
    return rtn;
}

char *region_to_string(const bed_region r)
{
    return fc_asprintf("region string", "%s:%" PRId64 "-%" PRId64, r->chr, r->start, r->end);
}

int region_from_bed(FILE *bed_fp, char **chr, int *start, int *end)
{
    char *line = NULL;
    char *line_copy = NULL;
    size_t len = 0;
    int rtn = 0;

    *start = -1;
    *end = -1;

    if (getline(&line, &len, bed_fp) != -1) {
        char *newline_pos = strchr(line, '\n');
        if (newline_pos != NULL) {
            *newline_pos = '\0'; // Null-terminate the string at the newline
        }
        line_copy = fc_strdup(line, "BED line for error reporting"); // Copy line for error reporting

        // get chromosome
        char *tok = strtok(line, "\t");
        if (tok == NULL) {
            FC_LOG_ERROR("Missing chromosome field in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
        size_t tok_len = strlen(tok);
        *chr = fc_realloc(*chr, (tok_len + 1) * sizeof(char), "chr");
        memcpy(*chr, tok, tok_len + 1);

        // get start coordinate
        tok = strtok(NULL, "\t");
        if (tok == NULL) {
            FC_LOG_ERROR("Missing start field in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
        char *endptr;
        errno = 0;
        long start_val = strtol(tok, &endptr, 10);
        if (*endptr != '\0' || errno == ERANGE || start_val < 0 || start_val > INT_MAX) {
            FC_LOG_ERROR("Invalid or out-of-range start coordinate in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
        *start = (int)start_val;

        // get end coordinate
        tok = strtok(NULL, "\t");
        if (tok == NULL) {
            FC_LOG_ERROR("Missing end field in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
        errno = 0;
        long end_val = strtol(tok, &endptr, 10);
        if (*endptr != '\0' || errno == ERANGE || end_val < 0 || end_val > INT_MAX) {
            FC_LOG_ERROR("Invalid or out-of-range end coordinate in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
        *end = (int)end_val;

        // Validate start and end
        if (*start < 0 || *end < 0 || *start >= *end) {
            FC_LOG_ERROR("Invalid region in BED file line: '%s'.", line_copy);
            rtn = -2;
            goto cleanup;
        }
    } else {
        rtn = -1; // EOF
    }

cleanup:
    free(line);
    free(line_copy);
    return rtn;
}

// Initialize the region iterator
regiter init_region_iterator(const char *bed_file, const char *single_region, const sam_hdr_t *hdr)
{
    regiter it = { 0 };
    it.hdr = hdr;
    if (bed_file != NULL) {
        it.bed_fp = fopen(bed_file, "r");
        if (it.bed_fp == NULL) {
            FC_LOG_ERROR("Unable to open BED file: %s", bed_file);
            it.error |= 1;
        }
    } else if (single_region != NULL) {
        it.single_region = fc_strdup(single_region, "single region string");
        it.mode = 1; // we'll process this first then switch to BED file mode
    }
    return it;
}

// Clean up the iterator
void destroy_region_iterator(regiter *it)
{
    if (it->chr != NULL) {
        free(it->chr);
    }
    if (it->single_region != NULL) {
        free(it->single_region);
    }
    if (it->bed_fp != NULL) {
        fclose(it->bed_fp);
    }
}

// Get the next region.
regiter_status next_region(regiter *it)
{
    if (it == NULL) {
        return REGITER_EOF;
    }
    regiter_status rtn = REGITER_OK;
    if (it->mode == 0) {
        if (it->bed_fp == NULL) {
            return REGITER_EOF;
        }
        rtn = (regiter_status)region_from_bed(it->bed_fp, &it->chr, &it->start, &it->end);
    } else if (it->mode == 1) {
        it->mode = (it->bed_fp != NULL) ? 0 : 2;
        if (it->single_region) {
            rtn = (regiter_status)region_from_string(it->single_region, &it->chr, &it->start, &it->end);
            if (rtn == REGITER_EOF) {
                FC_LOG_WARN("Failed to parse region string: %s", it->single_region);
                it->error |= 2;
                rtn = REGITER_PARSE_ERROR;
            }
        } else {
            rtn = REGITER_EOF;
        }
    } else {
        return REGITER_EOF;
    }

    if (rtn == REGITER_OK) {
        // check reference exists and tidy up length
        int tid = sam_hdr_name2tid((sam_hdr_t *)it->hdr, it->chr);
        if (tid < 0) {
            FC_LOG_WARN("Failed to find reference '%s' in BAM header.", it->chr);
            rtn = REGITER_REF_MISSING;
        } else {
            it->tid = tid;
            size_t ref_length = (size_t)sam_hdr_tid2len(it->hdr, tid);
            int ns = min(it->start, (int)ref_length);
            int ne = min(it->end, (int)ref_length);
            if (ns >= ne) {
                FC_LOG_WARN("Zero-length region created after truncating to reference length (%zu) '%s:%d-%d'.", ref_length, it->chr, it->start, it->end);
                rtn = REGITER_PARSE_ERROR;
            } else {
                it->start = ns;
                it->end = ne;
            }
            it->n_regions++;
        }
    }
    return rtn;
}

// Sorting regions
static int bed_region_cmp(const void *pa, const void *pb)
{
    const _bed_region *a = (const _bed_region *)pa;
    const _bed_region *b = (const _bed_region *)pb;

    if (a->tid != b->tid) {
        return (a->tid < b->tid) ? -1 : 1;
    }
    if (a->start != b->start) {
        return (a->start < b->start) ? -1 : 1;
    }
    if (a->end != b->end) {
        return (a->end < b->end) ? -1 : 1;
    }
    return 0;
}

void bed_regions_sort_by_header(bed_regions br)
{
    if (!br || br->n_regions <= 1) {
        return;
    }
    qsort(br->regions, br->n_regions,
        sizeof(br->regions[0]), bed_region_cmp);
}

// Parse a complete BED file into memory.
// Skips invalid/unknown regions (next_region already warns).
bed_regions init_bed(const char *bed_path, const sam_hdr_t *hdr)
{
    if (NULL == bed_path || NULL == hdr) {
        return NULL;
    }

    regiter it = init_region_iterator(bed_path, NULL, hdr);
    if (it.bed_fp == NULL) {
        exit(EXIT_FAILURE);
    }

    _bed_regions *out = fc_alloc(1, sizeof(_bed_regions), "bed_regions");
    out->_capacity = 1024;
    out->regions = fc_alloc(out->_capacity, sizeof(*out->regions), "bed_region[]");
    out->n_regions = 0;

    for (;;) {
        regiter_status rc = next_region(&it);
        if (rc == REGITER_OK) {
            if (out->n_regions == out->_capacity) {
                size_t new_cap = out->_capacity * 2;
                out->regions = fc_realloc(
                    out->regions, new_cap * sizeof(*out->regions), "bed_region[]");
                out->_capacity = new_cap;
            }
            bed_region dst = &out->regions[out->n_regions];
            dst->chr = fc_strdup(it.chr, "region chromosome name");
            dst->start = it.start;
            dst->end = it.end;
            dst->tid = it.tid;
            out->n_regions++;
            continue;
        }
        if (rc == REGITER_EOF) {
            break;
        }
        // Parse/ref errors are already warned by next_region(); skip.
    }

    destroy_region_iterator(&it);

    // shrink-to-fit
    if (out->n_regions < out->_capacity) {
        out->regions = fc_realloc(
            out->regions, out->n_regions * sizeof(*out->regions), "bed_region[]");
        out->_capacity = out->n_regions;
    }

    bed_regions_sort_by_header(out);
    return out;
}

bed_regions init_bed_from_sam(const sam_hdr_t *hdr, uint32_t segment_length)
{
    if (hdr == NULL) {
        return NULL;
    }

    // Create a bed_regions structure
    _bed_regions *regions = fc_alloc(1, sizeof(_bed_regions), "bed_regions");
    regions->_capacity = 1024;
    regions->regions = fc_alloc(regions->_capacity, sizeof(*regions->regions), "bed_region[]");
    regions->n_regions = 0;

    // Iterate through the header lines
    for (int i = 0; i < hdr->n_targets; ++i) {
        uint32_t span = segment_length == 0 ? hdr->target_len[i] : segment_length;

        for (uint32_t start = 0; start < hdr->target_len[i]; start += span) {
            if (regions->n_regions == regions->_capacity) {
                size_t new_cap = regions->_capacity * 2;
                regions->regions = fc_realloc(
                    regions->regions, new_cap * sizeof(*regions->regions), "bed_region[]");
                regions->_capacity = new_cap;
            }

            bed_region dst = &regions->regions[regions->n_regions];
            dst->chr = fc_strdup(hdr->target_name[i], "whole genome region chromosome");
            dst->start = start;
            dst->end = min(start + span, (uint32_t)hdr->target_len[i]);
            dst->tid = i;
            regions->n_regions++;
        }
    }

    // shrink-to-fit
    if (regions->n_regions < regions->_capacity) {
        regions->regions = (bed_region)fc_realloc(
            regions->regions, regions->n_regions * sizeof(*regions->regions), "bed_region[]");
        regions->_capacity = regions->n_regions;
    }

    return regions;
}

void destroy_bed(bed_regions regions)
{
    if (regions == NULL) {
        return;
    }
    if (NULL != regions->regions) {
        for (size_t i = 0; i < regions->n_regions; ++i) {
            free(regions->regions[i].chr);
        }
        free(regions->regions);
    }
    free(regions);
}
