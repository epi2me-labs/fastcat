#include "cli/fastlint/main.h"

#include <stdio.h>
#include <stdlib.h>

#include "args.h"
#include "core/fc_log.h"
#include "core/kseq_gz.h"
#include "vendor/sdust/sdust.h"

int fastlint_main(int argc, char *argv[])
{
    fc_fastq_stream *fp;
    kseq_t *ks;

    fc_log_config log_cfg = {
        .level = FC_LOG_INFO,
        .module = "FASTLINT",
        .target = NULL,
        .color = FC_LOG_COLOR_AUTO
    };
    fc_log_init(&log_cfg);

    arguments_t args = fastlint_parse_arguments(argc, argv);

    for (size_t i = 0; i < args.nfiles; ++i) {
        const char *input_name = args.fastq[i];
        if ((strcmp(args.fastq[i], "-") == 0) && (strlen(args.fastq[i]) == 1)) {
            fp = fc_fastq_stream_open_fd(fileno(stdin));
        } else {
            fp = fc_fastq_stream_open_path(args.fastq[i], 1);
        }
        if (fp == NULL) {
            FC_LOG_ERROR("Failed to open FASTQ input: %s", input_name);
            fc_log_destroy();
            return EXIT_FAILURE;
        }
        ks = kseq_init(fp);
        if (ks == NULL) {
            FC_LOG_ERROR("Failed to initialize FASTQ reader for input: %s", input_name);
            (void)fc_fastq_stream_close(fp);
            fc_log_destroy();
            return EXIT_FAILURE;
        }

        int kseq_status = 0;
        while ((kseq_status = kseq_read(ks)) >= 0) {
            uint64_t *r;
            int n;
            r = sdust(0, (uint8_t *)ks->seq.s, -1, args.t, args.w, &n);

            // Calculate total masked length
            int masked_bases = 0;
            for (int i = 0; i < n; ++i) {
                int start = (int)(r[i] >> 32);
                int end = (int)r[i];
                masked_bases += (end - start);
            }
            int read_length = ks->seq.l;
            double masked_fraction = read_length == 0 ? 0.0 : (double)masked_bases / read_length;
            if (masked_fraction > args.max_proportion) {
                FC_LOG_WARN("Read %s masked fraction %.2f exceeds threshold %.2f, skipping.",
                    ks->name.s, masked_fraction, args.max_proportion);
            } else {
                // for reasons unknown, we care about preserving the separator
                // between the name and the comment, but we don't have that. If
                // there are tabs in the comment, there's a good chance the
                // separator was a tab, particularly if data went through
                // fastcat, which emits SAM-style header tags.
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
        if (kseq_status < -1) {
            FC_LOG_ERROR("Failed to parse FASTQ input: %s", input_name);
            kseq_destroy(ks);
            (void)fc_fastq_stream_close(fp);
            fc_log_destroy();
            return EXIT_FAILURE;
        }
        kseq_destroy(ks);
        if (fc_fastq_stream_close(fp) != 0) {
            FC_LOG_ERROR("Failed to close FASTQ input stream: %s", input_name);
            fc_log_destroy();
            return EXIT_FAILURE;
        }
    }
    fc_log_destroy();
    return 0;
}
