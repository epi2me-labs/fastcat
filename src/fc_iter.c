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
#include <dirent.h>

#include "fc_log.h"
#include "fc_bam_utils.h"
#include "fc_filter.h"

void fc_file_stats_init(fc_file_stats *stats) {
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

static int fc_bam_iter_set_region(fc_bam_iter *iter, const char *chr, hts_pos_t start, hts_pos_t end);

int fc_fastq_iter_open(fc_fastq_iter *iter, const char *fname, const fc_iter_opts *opts) {
    if (!iter || !fname || fname[0] == '\0') return EINVAL;
    memset(iter, 0, sizeof(*iter));
    iter->fname = fname;
    fc_file_stats_init(&iter->stats);

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
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        return errno;
    }

    if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        return EISDIR;
    }

    iter->fp = gzopen(fname, "r");
    if (iter->fp == NULL) {
        FC_LOG_ERROR("could not open file %s: %s", fname, strerror(errno));
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
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
    const fc_iter_opts *opts) {
    if (!iter || !fp || !hdr) return -1;
    memset(iter, 0, sizeof(*iter));
    fc_file_stats_init(&iter->stats);

    // default: no length/quality filtering
    iter->opts.min_length = 0;
    iter->opts.max_length = (size_t)-1;
    iter->opts.min_qscore = 0.0f;

    iter->fp = fp;
    iter->hdr = hdr;
    iter->idx = idx;
    iter->own_idx = false;

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

    const char *bed_path = opts ? opts->bed : NULL;
    const char *region_string = opts ? opts->region : NULL;
    iter->use_regions = (bed_path != NULL) || (region_string != NULL);

    if (iter->use_regions) {
        if (iter->idx == NULL && fp->fn != NULL) {
            iter->idx = sam_index_load(fp, fp->fn);
            iter->own_idx = (iter->idx != NULL);
        }
        if (iter->idx == NULL) {
            FC_LOG_ERROR("Cannot find index for BAM when regions are requested.");
            return -1;
        }

        iter->b = bam_init1();
        if (iter->b == NULL) {
            if (iter->own_idx) hts_idx_destroy(iter->idx);
            iter->idx = NULL;
            return -1;
        }

        iter->reg = init_region_iterator(bed_path, region_string, hdr);

        int check = next_region(&iter->reg);
        if (check == -1) {
            iter->reg_active = false;
            return 0;
        }
        if (check == -2 && bed_path == NULL && region_string != NULL) {
            FC_LOG_ERROR("Failed to parse region string: %s", region_string);
            return -1;
        }
        while (check != 0) {
            check = next_region(&iter->reg);
            if (check == -1) break;
        }
        iter->reg_active = (check == 0);
        if (iter->reg_active) {
            if (fc_bam_iter_set_region(iter, iter->reg.chr, iter->reg.start, iter->reg.end) != 0) {
                iter->reg_active = false;
            }
        }
        return 0;
    }

    iter->bam = fc_bam_iter_data_create(fp, idx, hdr, NULL, 0, INT64_MAX);
    if (iter->bam == NULL) return -1;

    iter->b = bam_init1();
    if (iter->b == NULL) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
        return -1;
    }

    return 0;
}

static int fc_bam_iter_set_region(fc_bam_iter *iter, const char *chr, hts_pos_t start, hts_pos_t end) {
    if (iter->bam) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
    }
    iter->bam = fc_bam_iter_data_create(iter->fp, iter->idx, iter->hdr, chr, start, end);
    return iter->bam ? 0 : -1;
}


fc_read *fc_bam_iter_next(fc_bam_iter *iter, const sam_hdr_t *hdr) {
    if (!iter || !iter->b) return NULL;

    for (;;) {
        // If no active bam iterator, try to advance to next region (if any)
        if (!iter->bam) {
            if (!iter->use_regions || !iter->reg_active) {
                iter->status = -1;
                return NULL;
            }
        }

        iter->status = iter->bam ? fc_bam_read(iter->bam, iter->b) : -1;

        if (iter->status < 0) {
            if (iter->use_regions) {
                // advance to next region
                int check = next_region(&iter->reg);
                while (check > 0) { // skip invalid regions
                    check = next_region(&iter->reg);
                }
                if (check == 0) {
                    iter->reg_active = true;
                    // reset iterator for new region
                    // hdr and idx are already set
                    if (iter->bam) {
                        fc_bam_iter_data_destroy(iter->bam);
                        iter->bam = NULL;
                    }
                    if (fc_bam_iter_set_region(iter, iter->reg.chr, iter->reg.start, iter->reg.end) != 0) {
                        iter->reg_active = false;
                        continue;
                    }
                    iter->status = fc_bam_read(iter->bam, iter->b);
                    if (iter->status < 0) {
                        continue; // empty region, move on
                    }
                } else {
                    iter->reg_active = false;
                    iter->status = -1;
                    return NULL;
                }
            } else {
                return NULL;
            }
        }

        if (iter->status < 0) {
            return NULL;
        }

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
    if (iter->use_regions) {
        destroy_region_iterator(&iter->reg);
    }
    if (iter->own_idx && iter->idx) {
        hts_idx_destroy(iter->idx);
        iter->idx = NULL;
    }
}

// ---- Unified file iterator (FASTQ + BAM) ----

static bool has_suffix(const char *path, const char *suffix) {
    size_t lp = strlen(path);
    size_t ls = strlen(suffix);
    if (lp < ls) return false;
    return strcmp(path + (lp - ls), suffix) == 0;
}

static bool is_fastq_ext(const char *path) {
    return has_suffix(path, ".fastq") ||
           has_suffix(path, ".fastq.gz") ||
           has_suffix(path, ".fq") ||
           has_suffix(path, ".fq.gz");
}

static bool is_bam_ext(const char *path) {
    return has_suffix(path, ".bam") || has_suffix(path, ".sam");
}

static bool classify_kind(const char *path, fc_file_kind *out_kind) {
    if (is_bam_ext(path)) {
        if (out_kind) *out_kind = FC_FILE_BAM;
        return true;
    }
    if (is_fastq_ext(path)) {
        if (out_kind) *out_kind = FC_FILE_FASTQ;
        return true;
    }
    return false;
}

static int append_file(fc_iter *it, const char *path, fc_file_kind kind) {
    size_t new_count = it->n_files + 1;
    fc_file_entry *tmp = realloc(it->files, new_count * sizeof(*tmp));
    if (!tmp) return -1;
    it->files = tmp;
    it->files[it->n_files].path = strdup(path);
    if (!it->files[it->n_files].path) return -1;
    it->files[it->n_files].kind = kind;
    it->n_files = new_count;
    return 0;
}

typedef struct {
    char *path;
    int depth;
} pending_path;

int fc_iter_collect_paths(
    char **paths,
    size_t n_paths,
    int recurse,
    bool verbose,
    fc_file_entry **out_files,
    size_t *out_n_files) {
    if (!paths || n_paths == 0 || !out_files || !out_n_files) return EINVAL;
    (void)verbose; // logging occurs when iterating, not during collection

    int status = 0;
    fc_file_entry *files = NULL;
    size_t n_files = 0;

    pending_path *queue = NULL;
    size_t qlen = 0, qcap = 0, qhead = 0;

    // seed queue
    for (size_t i = 0; i < n_paths; ++i) {
        if (qlen == qcap) {
            size_t newcap = qcap ? qcap * 2 : 8;
            pending_path *tmp = realloc(queue, newcap * sizeof(*tmp));
            if (!tmp) { status = ENOMEM; goto done; }
            queue = tmp;
            qcap = newcap;
        }
        queue[qlen].path = strdup(paths[i]);
        queue[qlen].depth = recurse;
        if (!queue[qlen].path) { status = ENOMEM; goto done; }
        qlen++;
    }

    while (qhead < qlen) {
        pending_path cur = queue[qhead++];
        struct stat st = {0};
        if (stat(cur.path, &st) != 0) {
            FC_LOG_ERROR("could not stat '%s': %s", cur.path, strerror(errno));
            status = (status == 0) ? errno : status;
            free(cur.path);
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            if (cur.depth <= 0) {
                free(cur.path);
                continue;
            }
            DIR *dir = opendir(cur.path);
            if (!dir) {
                FC_LOG_ERROR("could not open directory '%s': %s", cur.path, strerror(errno));
                status = (status == 0) ? errno : status;
                free(cur.path);
                continue;
            }
            struct dirent *ent;
            while ((ent = readdir(dir)) != NULL) {
                if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) continue;
                size_t clen = strlen(cur.path) + 1 + strlen(ent->d_name) + 1;
                char *child = malloc(clen);
                if (!child) { status = (status == 0) ? ENOMEM : status; continue; }
                snprintf(child, clen, "%s/%s", cur.path, ent->d_name);
                if (qlen == qcap) {
                    size_t newcap = qcap ? qcap * 2 : 8;
                    pending_path *tmp = realloc(queue, newcap * sizeof(*tmp));
                    if (!tmp) { free(child); status = (status == 0) ? ENOMEM : status; continue; }
                    queue = tmp;
                    qcap = newcap;
                }
                queue[qlen].path = child;
                queue[qlen].depth = cur.depth - 1;
                qlen++;
            }
            closedir(dir);
            free(cur.path);
            continue;
        }

        fc_file_kind kind;
        if (!classify_kind(cur.path, &kind)) {
            free(cur.path);
            continue; // skip unrecognised file
        }
        fc_iter tmp = {.files = files, .n_files = n_files};
        if (append_file(&tmp, cur.path, kind) != 0) {
            status = (status == 0) ? ENOMEM : status;
            free(cur.path);
            continue;
        }
        files = tmp.files;
        n_files = tmp.n_files;
        free(cur.path);
    }

done:
    for (size_t i = qhead; i < qlen; ++i) {
        free(queue[i].path);
    }
    free(queue);
    if (status != 0) {
        for (size_t i = 0; i < n_files; ++i) free(files[i].path);
        free(files);
        files = NULL;
        n_files = 0;
    }
    *out_files = files;
    *out_n_files = n_files;
    return status;
}

int fc_iter_init(fc_iter *it, char **paths, size_t n_paths, int recurse, bool verbose) {
    if (!it || !paths || n_paths == 0) return EINVAL;
    memset(it, 0, sizeof(*it));
    it->recurse = recurse;
    it->verbose = verbose;

    return fc_iter_collect_paths(paths, n_paths, recurse, verbose, &it->files, &it->n_files);
}

bool fc_iter_next_file(fc_iter *it, const fc_file_entry **out_file) {
    if (!it || it->idx >= it->n_files) return false;
    it->current = &it->files[it->idx++];
    if (it->verbose) {
        FC_LOG_INFO("Processing %s", it->current->path);
    }
    if (out_file) *out_file = it->current;
    return true;
}

int fc_iter_open_record_cursor(
    fc_iter *it,
    const fc_iter_opts *opts) {
    if (!it || !it->current) return EINVAL;
    it->cursor_open = false;
    if (it->current->kind == FC_FILE_FASTQ) {
        fc_file_stats_init(&it->fqi.stats);
        int rc = fc_fastq_iter_open(&it->fqi, it->current->path, opts);
        if (rc != 0) {
            it->fqi.stats.failures[F_UNKNOWN_ERROR] = 1;
        }
        it->cursor_open = (rc == 0);
        return rc;
    }
    // BAM path
    fc_file_stats_init(&it->bmi.stats);
    it->bam_fp = hts_open(it->current->path, "rb");
    if (!it->bam_fp) {
        it->bmi.stats.failures[F_UNKNOWN_ERROR] = 1;
        return errno;
    }
    it->bam_hdr = sam_hdr_read(it->bam_fp);
    if (!it->bam_hdr) {
        hts_close(it->bam_fp);
        it->bam_fp = NULL;
        it->bmi.stats.failures[F_UNKNOWN_ERROR] = 1;
        return EINVAL;
    }
    int rc = fc_bam_iter_open(&it->bmi, it->bam_fp, NULL, it->bam_hdr, opts);
    if (rc != 0) {
        sam_hdr_destroy(it->bam_hdr);
        it->bam_hdr = NULL;
        hts_close(it->bam_fp);
        it->bam_fp = NULL;
        it->bmi.stats.failures[F_UNKNOWN_ERROR] = 1;
        return rc;
    }
    it->cursor_open = true;
    return 0;
}

fc_read *fc_iter_next_record(fc_iter *it, const sam_hdr_t *hdr) {
    if (!it || !it->current) return NULL;
    if (it->current->kind == FC_FILE_FASTQ) {
        return fc_fastq_iter_next(&it->fqi);
    }
    const sam_hdr_t *use_hdr = hdr ? hdr : it->bam_hdr;
    return fc_bam_iter_next(&it->bmi, use_hdr);
}

int fc_iter_finish_file(fc_iter *it) {
    if (!it || !it->current) return EINVAL;
    int rc = 0;
    if (it->current->kind == FC_FILE_FASTQ) {
        rc = fc_fastq_iter_finish(&it->fqi, NULL);
    } else {
        rc = fc_bam_iter_finish(&it->bmi, NULL);
        fc_bam_iter_close(&it->bmi);
        if (it->bam_hdr) {
            sam_hdr_destroy(it->bam_hdr);
            it->bam_hdr = NULL;
        }
        if (it->bam_fp) {
            hts_close(it->bam_fp);
            it->bam_fp = NULL;
        }
    }
    it->cursor_open = false;
    return rc;
}

void fc_iter_close(fc_iter *it) {
    if (!it) return;
    if (it->cursor_open && it->current) {
        if (it->current->kind == FC_FILE_BAM) {
            fc_bam_iter_close(&it->bmi);
            if (it->bam_hdr) sam_hdr_destroy(it->bam_hdr);
            if (it->bam_fp) hts_close(it->bam_fp);
        } else {
            fc_fastq_iter_finish(&it->fqi, NULL);
        }
    }
    for (size_t i = 0; i < it->n_files; ++i) {
        free(it->files[i].path);
    }
    free(it->files);
    it->files = NULL;
    it->n_files = 0;
    it->idx = 0;
    it->current = NULL;
}

void fc_iter_merge_file_stats(const fc_iter *it, fc_file_stats *dst) {
    if (!it || !dst) return;
    const fc_file_stats *src = fc_iter_stats(it);
    if (!src) return;
    double qsum = dst->mean_quality * (double)dst->n_seqs + src->mean_quality * (double)src->n_seqs;
    dst->n_seqs += src->n_seqs;
    dst->n_bases += src->n_bases;
    if (src->n_seqs > 0) {
        dst->min_length = (dst->min_length == 0 || dst->min_length == SIZE_MAX)
            ? src->min_length
            : min(dst->min_length, src->min_length);
        dst->max_length = max(dst->max_length, src->max_length);
    }
    for (size_t i = 0; i < NUM_FAILURE_CODES; ++i) {
        dst->failures[i] += src->failures[i];
    }
    if (dst->n_seqs > 0) {
        dst->mean_quality = qsum / (double)dst->n_seqs;
    }
}

const fc_file_stats *fc_iter_stats(const fc_iter *it) {
    if (!it || !it->current) return NULL;
    if (it->current->kind == FC_FILE_FASTQ) {
        return &it->fqi.stats;
    }
    return &it->bmi.stats;
}
