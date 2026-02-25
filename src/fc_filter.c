#include "fc_filter.h"

#include <errno.h>
#include <string.h>

#include "fc_bam_utils.h"
#include "fc_dust.h"

bool fc_filter_dust(const uint8_t *seq, size_t len, int t, int w, double max_fraction)
{
    return fc_dust_fraction(seq, len, t, w) <= max_fraction;
}

bool fc_filter_int_tag(const bam1_t *b, const char tag_name[2], int tag_value, bool keep_missing)
{
    if (!tag_name || tag_name[0] == '\0') {
        return true;
    }
    uint8_t *tag = fc_bam_get_tag_caseinsensitive(b, tag_name);
    if (tag == NULL) {
        return keep_missing;
    }
    errno = 0;
    int value = bam_aux2i(tag);
    if (value == 0 && errno == EINVAL) {
        return false;
    }
    return value == tag_value;
}

bool fc_filter_rg(const bam1_t *b, const char *read_group, bool keep_missing)
{
    if (read_group == NULL) {
        return true;
    }
    uint8_t *rg = fc_bam_get_tag_caseinsensitive(b, "RG");
    if (rg == NULL) {
        return keep_missing;
    }
    errno = 0;
    const char *value = bam_aux2Z(rg);
    if (value == NULL || (value[0] == '\0' && errno == EINVAL)) {
        return false;
    }
    return strcmp(read_group, value) == 0;
}
