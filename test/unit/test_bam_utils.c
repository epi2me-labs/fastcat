#include <string.h>
#include <stdlib.h>

#include <zlib.h>

#include "unity.h"
#include "htslib/kseq.h"

#ifndef KSEQ_DECLARED
#define KSEQ_DECLARED
KSEQ_DECLARE(gzFile)
#endif

#include "fastcat/bam_utils.h"
#include "fastqcomments.h"

static void test_kseq_to_bam_unmapped_aux(void) {
    kseq_t seq = {0};
    seq.name.s = strdup("read1");
    seq.name.l = strlen(seq.name.s);
    seq.seq.s = strdup("ACGT");
    seq.seq.l = strlen(seq.seq.s);
    seq.qual.s = strdup("####"); // ASCII 35 -> phred 2 after offset
    seq.qual.l = strlen(seq.qual.s);

    kstring_t comment = {0};
    read_meta meta = create_read_meta(&comment);
    ksprintf(meta->tags_str, "XY:Z:foo");

    bam1_t* bam = kseq_to_bam(&seq, meta);
    TEST_ASSERT_TRUE(bam->core.flag & BAM_FUNMAP);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.tid);
    TEST_ASSERT_EQUAL_INT(-1, bam->core.pos);
    TEST_ASSERT_EQUAL_INT(4, bam->core.l_qseq);
    TEST_ASSERT_EQUAL_STRING("read1", bam_get_qname(bam));

    uint8_t* qual = bam_get_qual(bam);
    TEST_ASSERT_EQUAL_UINT8(2, qual[0]); // 35 - 33

    uint8_t* aux = bam_aux_get(bam, "XY");
    TEST_ASSERT_NOT_NULL(aux);
    const char* aux_str = bam_aux2Z(aux);
    TEST_ASSERT_NOT_NULL(aux_str);
    TEST_ASSERT_EQUAL_STRING("foo", aux_str);

    bam_destroy1(bam);
    destroy_read_meta(meta);
    free(seq.name.s);
    free(seq.seq.s);
    free(seq.qual.s);
}

void bam_utils_suite(void) {
    RUN_TEST(test_kseq_to_bam_unmapped_aux);
}
