#include "core/fc_dust.h"

#include <stdlib.h>

#include "vendor/sdust/sdust.h"

double fc_dust_fraction(const uint8_t *seq, size_t len, int t, int w)
{
    if (len == 0) {
        return 0.0;
    }
    int n = 0;
    uint64_t *ranges = sdust(0, (uint8_t *)seq, -1, t, w, &n);
    int masked_bases = 0;
    for (int i = 0; i < n; ++i) {
        int start = (int)(ranges[i] >> 32);
        int end = (int)ranges[i];
        masked_bases += (end - start);
    }
    free(ranges);
    return (double)masked_bases / (double)len;
}
