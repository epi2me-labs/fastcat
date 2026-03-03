#include <zlib.h>

#include "htslib/kseq.h"

// Single global instantiation of kseq for gzFile/gzread. Keep this in one
// translation unit only.
KSEQ_INIT2(, gzFile, gzread)
