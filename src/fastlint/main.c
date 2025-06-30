
#include <zlib.h>
#include <stdio.h>
#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KSEQ_DECLARED

#include "../sdust/sdust.h"
#include "args.h"

int main(int argc, char *argv[]) {
    gzFile fp;
    kseq_t *ks;

    arguments_t args = parse_arguments(argc, argv);

    for (size_t i=0; i<args.nfiles; ++i) { 
        if ((strcmp(args.fastq[i], "-") == 0) && (strlen(args.fastq[i]) == 1)) {
            fp = gzdopen(fileno(stdin), "r");
        }
        else {
            fp = gzopen(args.fastq[i], "r");
        }
        ks = kseq_init(fp);

        while (kseq_read(ks) >= 0) {
            uint64_t *r;
            int n;
            r = sdust(0, (uint8_t*)ks->seq.s, -1, args.t, args.w, &n);
            
            // Calculate total masked length
            int masked_bases = 0;
            for (int i = 0; i < n; ++i) {
                int start = (int)(r[i] >> 32);
                int end = (int)r[i];
                masked_bases += (end - start);
            }
            int read_length = ks->seq.l;
            double masked_fraction = (double)masked_bases / read_length;
            if (masked_fraction > args.max_proportion) {
                fprintf(stderr, "Read %s masked fraction %.2f exceeds threshold %.2f, skipping.\n",
                    ks->name.s, masked_fraction, args.max_proportion);
            }
            else {
                // for reasons unknown, we care about preserving the separator
                // between the name and the comment, but we don't have that. If
                // there are tabs in the comment, there's a good chance the
                // separator was a tab, particular if data went through fastcat
                // with the --reheader option.
                if (ks->comment.l > 0) {
                    char sep = ' ';
                    if (strchr(ks->comment.s, '\t') != NULL) {
                        sep = '\t';
                    }
                    printf("@%s%c%s\n%s\n+\n%s\n", ks->name.s, sep, ks->comment.s, ks->seq.s, ks->qual.s);
                } else {
                    printf("@%s\n%s\n+\n%s\n", ks->name.s, ks->seq.s, ks->qual.s);
                }
            }
            free(r);
        }
        kseq_destroy(ks);
        gzclose(fp);
    }
	return 0;
}
