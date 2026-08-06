#include "core/fc_regiter.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "core/fc_common.h"
#include "core/fc_log.h"

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

    if (fc_getline(&line, &len, bed_fp) != -1) {
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

// Load either a BED file or a single region into a normalised region collection.
static bed_regions load_bed_regions(const char *bed_path, const char *single_region,
    const sam_hdr_t *hdr, int *error);
static regiter_status normalise_region(
    const sam_hdr_t *hdr, const char *chr, int *start, int *end, int *tid);

regiter init_region_iterator(const char *bed_file, const char *single_region, const sam_hdr_t *hdr)
{
    regiter it = { 0 };
    it.hdr = hdr;
    it.regions = load_bed_regions(bed_file, single_region, hdr, &it.error);
    if (it.regions == NULL && bed_file != NULL) {
        FC_LOG_ERROR("Unable to open BED file: %s", bed_file);
        it.error |= 1;
    }
    return it;
}

// Clean up the iterator
void destroy_region_iterator(regiter *it)
{
    destroy_bed(it->regions);
}

static regiter_status normalise_region(
    const sam_hdr_t *hdr, const char *chr, int *start, int *end, int *tid)
{
    *tid = sam_hdr_name2tid((sam_hdr_t *)hdr, chr);
    if (*tid < 0) {
        FC_LOG_WARN("Failed to find reference '%s' in BAM header.", chr);
        return REGITER_REF_MISSING;
    }

    size_t ref_length = (size_t)sam_hdr_tid2len(hdr, *tid);
    int normalised_start = min(*start, (int)ref_length);
    int normalised_end = min(*end, (int)ref_length);
    if (normalised_start >= normalised_end) {
        FC_LOG_WARN("Zero-length region created after truncating to reference length (%zu) '%s:%d-%d'.", ref_length, chr, *start, *end);
        return REGITER_PARSE_ERROR;
    }
    *start = normalised_start;
    *end = normalised_end;
    return REGITER_OK;
}

// Get the next region.
regiter_status next_region(regiter *it)
{
    if (it == NULL || it->regions == NULL || it->region_index >= it->regions->n_regions) {
        return REGITER_EOF;
    }
    const bed_region region = &it->regions->regions[it->region_index++];
    it->chr = region->chr;
    it->start = region->start;
    it->end = region->end;
    it->tid = region->tid;
    return REGITER_OK;
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

void bed_regions_merge_overlaps(bed_regions regions)
{
    if (regions == NULL || regions->n_regions < 2) {
        return;
    }

    size_t out = 1;
    for (size_t i = 1; i < regions->n_regions; ++i) {
        bed_region previous = &regions->regions[out - 1];
        bed_region current = &regions->regions[i];
        if (previous->tid == current->tid && current->start <= previous->end) {
            if (current->end > previous->end) {
                previous->end = current->end;
            }
            free(current->chr);
            continue;
        }
        if (out != i) {
            regions->regions[out] = regions->regions[i];
        }
        out++;
    }
    regions->n_regions = out;
}

static void bed_regions_append(
    bed_regions regions, char *chr, int start, int end, int tid)
{
    if (regions->n_regions == regions->_capacity) {
        size_t new_cap = regions->_capacity * 2;
        regions->regions = fc_realloc(
            regions->regions, new_cap * sizeof(*regions->regions), "bed_region[]");
        regions->_capacity = new_cap;
    }
    bed_region dst = &regions->regions[regions->n_regions++];
    dst->chr = chr;
    dst->start = start;
    dst->end = end;
    dst->tid = tid;
}

// Parse, validate, and normalise a BED file or one region into memory.
// Invalid or unknown regions are warned about and skipped.
static bed_regions load_bed_regions(const char *bed_path, const char *single_region,
    const sam_hdr_t *hdr, int *error)
{
    if (hdr == NULL || (bed_path == NULL && single_region == NULL)) {
        return NULL;
    }

    FILE *bed_fp = NULL;
    if (bed_path != NULL) {
        bed_fp = fopen(bed_path, "r");
        if (bed_fp == NULL) {
            return NULL;
        }
    }

    _bed_regions *out = fc_alloc(1, sizeof(_bed_regions), "bed_regions");
    out->_capacity = 1024;
    out->regions = fc_alloc(out->_capacity, sizeof(*out->regions), "bed_region[]");
    out->n_regions = 0;

    if (bed_fp == NULL) {
        char *chr = NULL;
        int start = 0;
        int end = 0;
        int tid = -1;
        if (region_from_string((char *)single_region, &chr, &start, &end) != 0) {
            FC_LOG_WARN("Failed to parse region string: %s", single_region);
            if (error != NULL) {
                *error |= 2;
            }
        } else {
            regiter_status status = normalise_region(hdr, chr, &start, &end, &tid);
            if (status == REGITER_OK) {
                bed_regions_append(out, chr, start, end, tid);
                chr = NULL;
            }
            if (status == REGITER_PARSE_ERROR && error != NULL) {
                *error |= 2;
            }
        }
        free(chr);
    } else {
        for (;;) {
            char *chr = NULL;
            int start = 0;
            int end = 0;
            regiter_status rc = (regiter_status)region_from_bed(bed_fp, &chr, &start, &end);
            if (rc == REGITER_EOF) {
                free(chr);
                break;
            }
            if (rc == REGITER_OK) {
                int tid = -1;
                rc = normalise_region(hdr, chr, &start, &end, &tid);
                if (rc == REGITER_OK) {
                    bed_regions_append(out, chr, start, end, tid);
                    chr = NULL;
                }
            }
            free(chr);
            // Parse/ref errors are already warned; skip this BED line.
        }
    }
    if (bed_fp != NULL) {
        fclose(bed_fp);
    }

    // shrink-to-fit. realloc(ptr, 0) may return NULL after freeing ptr.
    if (out->n_regions == 0) {
        free(out->regions);
        out->regions = NULL;
        out->_capacity = 0;
    } else if (out->n_regions < out->_capacity) {
        out->regions = fc_realloc(
            out->regions, out->n_regions * sizeof(*out->regions), "bed_region[]");
        out->_capacity = out->n_regions;
    }

    bed_regions_sort_by_header(out);
    bed_regions_merge_overlaps(out);
    return out;
}

bed_regions init_bed(const char *bed_path, const sam_hdr_t *hdr)
{
    bed_regions out = load_bed_regions(bed_path, NULL, hdr, NULL);
    if (out == NULL) {
        FC_LOG_FATAL("Unable to open BED file: %s", bed_path);
    }
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
