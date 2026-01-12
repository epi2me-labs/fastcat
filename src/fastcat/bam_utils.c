#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/sam.h"

#include "bam_utils.h"

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


bam1_t* kseq_to_bam(kseq_t* seq, read_meta meta) {
    bam1_t* b = bam_init1();
    bam_set1(
        b,
        seq->name.l, seq->name.s,
        BAM_FUNMAP, -1, -1, 0,
        0, NULL,
        -1, -1, 0,
        seq->seq.l, seq->seq.s, seq->qual.s,
        0
    );
    // bam_set1 uses raw numeric qualities; convert from ASCII
    uint8_t* qual = bam_get_qual(b);
    for (int i = 0; i < b->core.l_qseq; ++i) {
        qual[i] -= 33;
    }

    if (meta != NULL && meta->tags_str != NULL && meta->tags_str->l > 0) {
        if (parse_and_set_aux_tags(b, meta->tags_str->s) < 0) {
            fprintf(stderr, "Error parsing auxiliary tags for read '%s'\n", seq->name.s);
            exit(1);
        }
    }

    return b;
}
