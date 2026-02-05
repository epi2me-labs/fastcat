#include <zlib.h>

#include "htslib/kseq.h"

KSEQ_INIT2(, gzFile, gzread)
#define KSEQ_DECLARED

#include "fc_iter.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "fc_log.h"
#include "fc_bam_utils.h"
#include "fc_filter.h"

static void init_file_stats(fc_file_stats *stats) {
    if (!stats) return;
    memset(stats, 0, sizeof(*stats));
    stats->min_length = SIZE_MAX;
}

static void finalize_file_stats(fc_file_stats *stats, double q_sum) {
    if (!stats) return;
    if (stats->n_seqs == 0) {
        stats->min_length = 0;
        stats->max_length = 0;
        stats->mean_quality = 0.0;
    } else {
        stats->mean_quality = q_sum / (double)stats->n_seqs;
        if (stats->min_length == SIZE_MAX) {
            stats->min_length = 0;
        }
    }
}

static inline void fc_iter_accum_add(size_t len, float mean_q,
        fc_file_stats *stats, double *q_sum, double *q_comp) {
    stats->n_seqs += 1;
    stats->n_bases += len;
    stats->min_length = (stats->min_length < len) ? stats->min_length : len;
    stats->max_length = (stats->max_length > len) ? stats->max_length : len;
    kahan_sum(q_sum, mean_q, q_comp);
}

static inline bool fc_iter_filter_length(size_t len, size_t min_len, size_t max_len, fc_file_stats *stats) {
    if (max_len != (size_t)-1 && len > max_len) {
        stats->failures[R_TOO_LONG]++;
        return false;
    }
    if (len < min_len) {
        stats->failures[R_TOO_SHORT]++;
        return false;
    }
    return true;
}

static inline bool fc_iter_filter_quality(float mean_q, float min_qscore, fc_file_stats *stats) {
    if (min_qscore <= 0.0f) return true;
    if (mean_q < min_qscore) {
        stats->failures[R_LOW_QUALITY]++;
        return false;
    }
    return true;
}

static float bam_iter_mean_quality(const bam1_t *b, bool force_recalc_qual) {
    if (!force_recalc_qual) {
        uint8_t *qs_tag = bam_aux_get(b, "QS");
        if (qs_tag != NULL) {
            errno = 0;
            float qs = bam_aux2f(qs_tag);
            if (!(qs == 0.0f && errno == EINVAL)) {
                return qs;
            }
        }
    }
    return mean_qual_from_bam_naive(bam_get_qual(b), b->core.l_qseq);
}

int fc_fastq_iter_open(fc_fastq_iter *iter, const char *fname, const fc_fastq_iter_opts *opts) {
    if (!iter || !fname || fname[0] == '\0') return EINVAL;
    memset(iter, 0, sizeof(*iter));
    iter->fname = fname;
    init_file_stats(&iter->stats);

    // default: no length/quality filtering
    iter->opts.min_length = 0;
    iter->opts.max_length = (size_t)-1;
    iter->opts.min_qscore = 0.0f;

    if (opts) {
        iter->opts = *opts;
        if (iter->opts.max_length == 0) iter->opts.max_length = (size_t)-1;
    }

    struct stat finfo;
    if (stat(fname, &finfo) == -1) {
        FC_LOG_ERROR("could not process file %s: %s", fname, strerror(errno));
        return errno;
    }

    if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
        return EISDIR;
    }

    iter->fp = gzopen(fname, "r");
    if (iter->fp == NULL) {
        FC_LOG_ERROR("could not open file %s: %s", fname, strerror(errno));
        return errno;
    }

    iter->seq = kseq_init(iter->fp);
    return 0;
}

fc_read *fc_fastq_iter_next(fc_fastq_iter *iter) {
    if (!iter || !iter->seq) return NULL;

    while ((iter->status = kseq_read(iter->seq)) != -1) {
        if (iter->status == -2) {
            iter->stats.failures[F_QUAL_TRUNCATED]++;
            iter->truncated = true;
            continue;
        } else if (iter->status == -3) {
            iter->stats.failures[F_STREAM_ERROR]++;
            return NULL;
        } else if (iter->status < 0) {
            iter->stats.failures[F_UNKNOWN_ERROR]++;
            return NULL;
        }
        if (iter->seq->qual.l == 0) {
            iter->stats.failures[F_QUAL_MISSING]++;
            iter->truncated = true;
            iter->status = -2;
            continue;
        }

        iter->truncated = false;
        iter->stats.failures[R_RECORD_OK]++;

        if (!fc_iter_filter_length(iter->seq->seq.l, iter->opts.min_length, iter->opts.max_length, &iter->stats)) {
            continue;
        }

        if (iter->opts.dust) {
            if (!fc_filter_dust((uint8_t *)iter->seq->seq.s, iter->seq->seq.l,
                    iter->opts.dust_t, iter->opts.dust_w, iter->opts.max_dust)) {
                iter->stats.failures[R_DUST_MASKED]++;
                continue;
            }
        }

        fc_read *fr = fc_read_from_kseq(iter->seq);
        if (!fc_iter_filter_quality(fr->meta->mean_q, iter->opts.min_qscore, &iter->stats)) {
            fc_read_destroy(fr);
            continue;
        }
        if (iter->opts.tag_name[0] != '\0') {
            if (!fc_filter_int_tag(fr->bam, iter->opts.tag_name, iter->opts.tag_value, false)) {
                fc_read_destroy(fr);
                continue;
            }
        }
        if (iter->opts.read_group != NULL) {
            if (!fc_filter_rg(fr->bam, iter->opts.read_group, false)) {
                fc_read_destroy(fr);
                continue;
            }
        }

        fc_iter_accum_add(iter->seq->seq.l, fr->meta->mean_q, &iter->stats, &iter->q_sum, &iter->q_comp);

        return fr;
    }

    if (iter->truncated) {
        iter->stats.failures[F_STREAM_ERROR]++;
        iter->status = -3;
    } else {
        iter->stats.failures[F_FILE_OK]++;
    }

    return NULL;
}

int fc_fastq_iter_finish(fc_fastq_iter *iter, fc_file_stats *out_stats) {
    if (!iter) return EINVAL;

    if (iter->stats.failures[F_STREAM_ERROR] > 0) {
        FC_LOG_WARN("file '%s' is possibly truncated.", iter->fname);
    } else if (iter->stats.failures[F_QUAL_MISSING] > 0) {
        FC_LOG_WARN("no quality string found for %" PRIu64 " records in file '%s'.",
            iter->stats.failures[F_QUAL_MISSING], iter->fname);
    } else if (iter->stats.failures[F_QUAL_TRUNCATED] > 0) {
        FC_LOG_WARN("truncated quality string found for %" PRIu64 " records in file '%s'.",
            iter->stats.failures[F_QUAL_TRUNCATED], iter->fname);
    } else if (iter->stats.failures[F_UNKNOWN_ERROR] > 0) {
        FC_LOG_WARN("unknown error reading file '%s'.", iter->fname);
    }

    finalize_file_stats(&iter->stats, iter->q_sum);

    if (out_stats) {
        *out_stats = iter->stats;
    }

    if (iter->seq) {
        kseq_destroy(iter->seq);
        iter->seq = NULL;
    }
    if (iter->fp) {
        gzclose(iter->fp);
        iter->fp = NULL;
    }

    return (iter->status == -1 || iter->status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

int fc_bam_iter_open(
    fc_bam_iter *iter,
    htsFile *fp,
    hts_idx_t *idx,
    sam_hdr_t *hdr,
    const char *chr,
    hts_pos_t start,
    hts_pos_t end,
    const fc_bam_iter_opts *opts) {
    if (!iter || !fp || !hdr) return -1;
    memset(iter, 0, sizeof(*iter));
    init_file_stats(&iter->stats);

    // default: no length/quality filtering
    iter->opts.min_length = 0;
    iter->opts.max_length = (size_t)-1;
    iter->opts.min_qscore = 0.0f;

    iter->bam = fc_bam_iter_data_create(fp, idx, hdr, chr, start, end);
    if (iter->bam == NULL) return -1;

    iter->b = bam_init1();
    if (iter->b == NULL) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
        return -1;
    }

    if (opts) {
        iter->opts = *opts;
        iter->has_opts = true;
        if (iter->opts.max_length == 0) iter->opts.max_length = (size_t)-1;
    }

    iter->apply_filters =
        iter->has_opts &&
        (iter->opts.min_length > 0 ||
         iter->opts.max_length != (size_t)-1 ||
         iter->opts.min_qscore > 0.0f ||
         iter->opts.dust ||
         iter->opts.tag_name[0] != '\0' ||
         iter->opts.read_group != NULL);

    return 0;
}

fc_read *fc_bam_iter_next(fc_bam_iter *iter, const sam_hdr_t *hdr) {
    if (!iter || !iter->bam || !iter->b) return NULL;

    while ((iter->status = fc_bam_read(iter->bam, iter->b)) >= 0) {
        iter->stats.failures[R_RECORD_OK]++;
        size_t len = iter->b->core.l_qseq;
        float mean_q = 0.0f;
        bool mean_q_set = false;

        if (iter->apply_filters) {
            if (!fc_iter_filter_length(len, iter->opts.min_length, iter->opts.max_length, &iter->stats)) {
                continue;
            }
            if (iter->opts.min_qscore > 0.0f) {
                mean_q = bam_iter_mean_quality(iter->b, iter->opts.force_recalc_qual);
                mean_q_set = true;
                if (!fc_iter_filter_quality(mean_q, iter->opts.min_qscore, &iter->stats)) {
                    continue;
                }
            }
            if (iter->opts.dust) {
                uint8_t *seq_data = bam_get_seq(iter->b);
                char *seq = xalloc(len + 1, sizeof(char), "bam dust seq");
                for (size_t i = 0; i < len; ++i) {
                    seq[i] = seq_nt16_str[bam_seqi(seq_data, i)];
                }
                seq[len] = '\0';
                bool pass = fc_filter_dust((uint8_t *)seq, len, iter->opts.dust_t,
                    iter->opts.dust_w, iter->opts.max_dust);
                free(seq);
                if (!pass) {
                    iter->stats.failures[R_DUST_MASKED]++;
                    continue;
                }
            }
            if (!fc_filter_int_tag(iter->b, iter->opts.tag_name, iter->opts.tag_value, false)) {
                continue;
            }
            if (!fc_filter_rg(iter->b, iter->opts.read_group, false)) {
                continue;
            }
        }

        if (!mean_q_set) {
            mean_q = bam_iter_mean_quality(iter->b, iter->opts.force_recalc_qual);
        }

        fc_iter_accum_add(len, mean_q, &iter->stats, &iter->q_sum, &iter->q_comp);

        return fc_read_from_bam(iter->b, hdr, iter->has_opts ? iter->opts.force_recalc_qual : false);
    }

    return NULL;
}

int fc_bam_iter_finish(fc_bam_iter *iter, fc_file_stats *out_stats) {
    if (!iter) return EINVAL;

    if (iter->status < -1) {
        iter->stats.failures[F_STREAM_ERROR]++;
    } else if (iter->status == -1) {
        iter->stats.failures[F_FILE_OK]++;
    } else if (iter->status < 0) {
        iter->stats.failures[F_UNKNOWN_ERROR]++;
    }

    finalize_file_stats(&iter->stats, iter->q_sum);

    if (out_stats) {
        *out_stats = iter->stats;
    }

    return (iter->status == -1 || iter->status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

void fc_bam_iter_close(fc_bam_iter *iter) {
    if (!iter) return;
    if (iter->bam) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
    }
    if (iter->b) {
        bam_destroy1(iter->b);
        iter->b = NULL;
    }
}
