#ifndef _BAM_INDEX_INDEX_H
#define _BAM_INDEX_INDEX_H


typedef struct bc_rec_t {
    size_t file_offset;
    size_t lqname;
    char *qname; 
} bc_rec_t;


typedef struct bc_idx_t {
    size_t version;
    size_t chunk_size;
    size_t n_chunks;
    size_t stored;  // TODO: disentangle header from contents
    bc_rec_t *recs; 
} bc_idx_t;


char* generate_index_filename(const char* input_bam, const char* input_index);
bc_idx_t *bc_idx_init();
bc_idx_t *bc_idx_init1();
void bc_idx_destroy(bc_idx_t *h);
bc_idx_t *bc_idx_read(FILE *fp);
int bc_idx_write_header(FILE* fp, bc_idx_t* idx);
int bc_idx_write(FILE* fp, bc_idx_t* idx, size_t offset, char* qname);

#endif