#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "unity.h"
#include "regiter.h"
#include "common.h"

static sam_hdr_t* make_header(void) {
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "100", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr2", "LN", "50", NULL);
    return hdr;
}

static void test_region_from_string_valid(void) {
    char* chr = NULL;
    int start = 0, end = 0;
    char* input = strdup("chr1:10-20");
    TEST_ASSERT_EQUAL_INT(0, region_from_string(input, &chr, &start, &end));
    TEST_ASSERT_EQUAL_STRING("chr1", chr);
    TEST_ASSERT_EQUAL_INT(9, start);   // hts_parse_reg converts to 0-based
    TEST_ASSERT_EQUAL_INT(20, end);    // half-open end
    free(chr);
    free(input);
}

static void test_region_from_string_error_return(void) {
    char* chr = NULL;
    int start = 0, end = 0;
    char* input = strdup("chr1:20-10"); // start >= end invalid
    TEST_ASSERT_EQUAL_INT(-1, region_from_string(input, &chr, &start, &end));
    TEST_ASSERT_NULL(chr);
    free(input);
}

static void test_region_from_string_invalid(void) {
    char* chr = NULL;
    int start = 0, end = 0;
    char* input = strdup("not_a_region");
    TEST_ASSERT_EQUAL_INT(0, region_from_string(input, &chr, &start, &end));
    TEST_ASSERT_EQUAL_STRING("not_a_region", chr);
    TEST_ASSERT_TRUE(start >= 0);
    TEST_ASSERT_TRUE(end >= start);
    free(chr);
    free(input);
}

static void test_region_to_string_roundtrip(void) {
    _bed_region region = {.chr = "chr2", .tid = 0, .start = 5, .end = 25};
    char* s = region_to_string(&region);
    TEST_ASSERT_EQUAL_STRING("chr2:5-25", s);

    char* chr = NULL;
    int start = 0, end = 0;
    TEST_ASSERT_EQUAL_INT(0, region_from_string(s, &chr, &start, &end));
    TEST_ASSERT_EQUAL_STRING("chr2", chr);
    TEST_ASSERT_EQUAL_INT(4, start); // 1-based to 0-based conversion inside region_from_string
    TEST_ASSERT_EQUAL_INT(25, end);

    free(chr);
    free(s);
}

static void test_region_from_bed_errors(void) {
    FILE* fp = tmpfile();
    fprintf(fp, "chr1\t5\n");              // missing end
    fprintf(fp, "chr1\tA\t10\n");         // bad start
    fprintf(fp, "chr1\t5\tX\n");          // bad end
    fprintf(fp, "chr1\t10\t5\n");         // start>=end
    rewind(fp);
    char* chr = NULL;
    int start = 0, end = 0;

    TEST_ASSERT_EQUAL_INT(-2, region_from_bed(fp, &chr, &start, &end)); // missing end
    TEST_ASSERT_EQUAL_INT(-2, region_from_bed(fp, &chr, &start, &end)); // bad start
    TEST_ASSERT_EQUAL_INT(-2, region_from_bed(fp, &chr, &start, &end)); // bad end
    TEST_ASSERT_EQUAL_INT(-2, region_from_bed(fp, &chr, &start, &end)); // start>=end
    TEST_ASSERT_EQUAL_INT(-1, region_from_bed(fp, &chr, &start, &end)); // EOF

    free(chr);
    fclose(fp);
}

static void test_region_from_bed_valid(void) {
    FILE* fp = tmpfile();
    fprintf(fp, "chr1\t0\t10\n");
    rewind(fp);

    char* chr = NULL;
    int start = -1, end = -1;
    TEST_ASSERT_EQUAL_INT(0, region_from_bed(fp, &chr, &start, &end));
    TEST_ASSERT_EQUAL_STRING("chr1", chr);
    TEST_ASSERT_EQUAL_INT(0, start);
    TEST_ASSERT_EQUAL_INT(10, end);

    free(chr);
    fclose(fp);
}

static void test_next_region_and_truncation(void) {
    sam_hdr_t* hdr = make_header();
    regiter it = init_region_iterator(NULL, "chr1:90-200", hdr);
    TEST_ASSERT_EQUAL_INT(0, next_region(&it));
    TEST_ASSERT_EQUAL_STRING("chr1", it.chr);
    TEST_ASSERT_EQUAL_INT(89, it.start); // truncated to ref length 100
    TEST_ASSERT_EQUAL_INT(100, it.end);
    destroy_region_iterator(&it);
    sam_hdr_destroy(hdr);
}

static void test_init_bed_from_sam_segments(void) {
    sam_hdr_t* hdr = make_header();
    bed_regions br = init_bed_from_sam(hdr, 50);
    TEST_ASSERT_EQUAL_size_t(3, br->n_regions); // chr1 split into two, chr2 one

    TEST_ASSERT_EQUAL_STRING("chr1", br->regions[0].chr);
    TEST_ASSERT_EQUAL_INT(0, br->regions[0].start);
    TEST_ASSERT_EQUAL_INT(50, br->regions[0].end);

    TEST_ASSERT_EQUAL_STRING("chr1", br->regions[1].chr);
    TEST_ASSERT_EQUAL_INT(50, br->regions[1].start);
    TEST_ASSERT_EQUAL_INT(100, br->regions[1].end);

    TEST_ASSERT_EQUAL_STRING("chr2", br->regions[2].chr);
    TEST_ASSERT_EQUAL_INT(0, br->regions[2].start);
    TEST_ASSERT_EQUAL_INT(50, br->regions[2].end);

    destroy_bed(br);
    sam_hdr_destroy(hdr);
}

void regiter_suite(void) {
    RUN_TEST(test_region_from_string_valid);
    RUN_TEST(test_region_from_string_error_return);
    RUN_TEST(test_region_from_string_invalid);
    RUN_TEST(test_region_to_string_roundtrip);
    RUN_TEST(test_region_from_bed_errors);
    RUN_TEST(test_region_from_bed_valid);
    RUN_TEST(test_next_region_and_truncation);
    RUN_TEST(test_init_bed_from_sam_segments);
}
