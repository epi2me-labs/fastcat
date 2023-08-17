#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "htslib/sam.h"

#include "index.h"

const size_t MAGIC_LEN = 5;
const char* FILE_MAGIC = "FANZ\0";


char* generate_index_filename(const char* input_bam, const char* input_index) {
    char* out_fn;

    if(input_index != NULL) {
        out_fn = malloc(strlen(input_index));
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_index);
    } else {
        out_fn = calloc(sizeof(char), strlen(input_bam) + 5);
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_bam);
        strcat(out_fn, ".bci");
    }
    return out_fn;
}

bc_idx_t *bc_idx_init(void) {
    bc_idx_t *h = (bc_idx_t*)calloc(1, sizeof(bc_idx_t));
    if (h == NULL) return NULL;
    // any init?
    return h;
}

bc_idx_t *bc_idx_init1(const size_t chunk_size) {
    bc_idx_t *idx = bc_idx_init();
    if (idx == NULL) return NULL;
    idx->version = 1;
    idx->chunk_size = chunk_size;
    idx->n_chunks = 0;
    idx->stored = 0;
    return idx;
}

void bc_idx_destroy(bc_idx_t *h) {
    if (h->stored > 0) {
        for (size_t i=0; i<h->stored; ++i) {
            if (h->recs[i].qname != NULL) {
                free(h->recs[i].qname);
            }
        }
        if (h->recs != NULL) {
            free(h->recs);
        }
    }
    free(h);
}

bc_idx_t *bc_idx_read(FILE *fp) {
    char buf[MAGIC_LEN];
    size_t magic_len = fread(&(buf), sizeof(char), MAGIC_LEN, fp);
    if (magic_len != MAGIC_LEN || memcmp(buf, FILE_MAGIC, MAGIC_LEN)) {
        fprintf(stderr, "Invalid BAM chunk index binary header.\n");
        return NULL;
    }
    bc_idx_t *h = bc_idx_init();
    if(h == NULL) {
        fprintf(stderr, "Failed to allocate header.\n");
        return NULL;
    }
    size_t items = 0;
    items += fread(&(h->version), sizeof(h->version), 1, fp);
    items += fread(&(h->chunk_size), sizeof(h->chunk_size), 1, fp);
    items += fread(&(h->n_chunks), sizeof(h->n_chunks), 1, fp);
    if (items != 3) {
        bc_idx_destroy(h);
        fprintf(stderr, "Invalid BAM chunk index binary header.\n");
        return NULL;
    }

    h->stored = h->n_chunks;
    h->recs = (bc_rec_t*)calloc(h->stored, sizeof(bc_rec_t));

    size_t valid = 0;
    char *msg = "Failed to read index contents. File is currupt.\n";
    for (size_t i=0; i<h->n_chunks; ++i, ++valid) {
        bc_rec_t *r = &(h->recs[i]);
        if (fread(&(r->file_offset), sizeof(r->file_offset), 1, fp) != 1) {
            fputs(msg, stderr); break;
        }
        if (fread(&(r->lqname), sizeof(r->lqname), 1, fp) != 1) {
            fputs(msg, stderr); break;
        }
        r->qname = (char*)calloc(r->lqname, sizeof(char));
        if (fread(r->qname, sizeof(char), r->lqname, fp) != r->lqname) {
            fputs(msg, stderr); break;
        }
    }
    if (valid != h->stored) {
        bc_idx_destroy(h);
        return NULL;
    }
    return h;
}

int bc_idx_write_header(FILE* fp, bc_idx_t* idx) {
    int rtn = 0;
    fseek(fp, 0, SEEK_SET);
    if (fwrite(FILE_MAGIC, sizeof(char), MAGIC_LEN, fp) != MAGIC_LEN) rtn = 1;
    if (fwrite(&(idx->version), sizeof(idx->version), 1, fp) != 1) rtn = 2;
    if (fwrite(&(idx->chunk_size), sizeof(idx->chunk_size), 1, fp) != 1) rtn = 3;
    if (fwrite(&(idx->n_chunks), sizeof(idx->n_chunks), 1, fp) != 1) rtn = 4;
    fseek(fp, 0, SEEK_END);
    return rtn;
}

int bc_idx_write(FILE* fp, bc_idx_t* idx, size_t offset, char* qname) {
    // write: file offset, length qname, qname
    size_t l_qname = strlen(qname) + 1;
    if (fwrite(&offset, sizeof(offset), 1, fp) != 1) return -1;
    if (fwrite(&l_qname, sizeof(l_qname), 1, fp) != 1) return -1; 
    if (fwrite(qname, sizeof(char), l_qname, fp) != l_qname) return -1;
    (idx->n_chunks)++;
    return (int)(idx->n_chunks);
}
