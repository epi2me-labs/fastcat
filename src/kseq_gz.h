#ifndef FASTCAT_KSEQ_GZ_H
#define FASTCAT_KSEQ_GZ_H

#include <zlib.h>

#include "htslib/kseq.h"

#ifndef KSEQ_DECLARED
KSEQ_DECLARE(gzFile)
#define KSEQ_DECLARED
#endif

#endif
