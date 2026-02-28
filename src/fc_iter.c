
#include "fc_iter.h"

#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fc_bam_utils.h"
#include "fc_filter.h"
#include "fc_input_list.h"
#include "fc_log.h"
#include "kseq_gz.h"

void fc_file_stats_init(fc_file_stats *stats)
{
    if (!stats) {
        return;
    }
    memset(stats, 0, sizeof(*stats));
    stats->min_length = SIZE_MAX;
}

static void finalize_file_stats(fc_file_stats *stats, double q_sum)
{
    if (!stats) {
        return;
    }
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

static inline void free_dust_buffer(dust_buffer *dust)
{
    if (dust->buf) {
        free(dust->buf);
        dust->buf = NULL;
        dust->size = 0;
    }
}

static inline void fc_iter_accum_add(size_t len, float mean_q,
    fc_file_stats *stats, double *q_sum, double *q_comp)
{
    stats->n_seqs += 1;
    stats->n_bases += len;
    stats->min_length = (stats->min_length < len) ? stats->min_length : len;
    stats->max_length = (stats->max_length > len) ? stats->max_length : len;
    kahan_sum(q_sum, mean_q, q_comp);
}

static inline bool fc_iter_filter_length(size_t len, size_t min_len, size_t max_len, fc_file_stats *stats)
{
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

static inline bool fc_iter_filter_quality(float mean_q, float min_qscore, fc_file_stats *stats)
{
    if (min_qscore <= 0.0f) {
        return true;
    }
    if (mean_q < min_qscore) {
        stats->failures[R_LOW_QUALITY]++;
        return false;
    }
    return true;
}

static void fc_iter_opts_init(fc_iter_opts *dst, const fc_iter_opts *src)
{
    if (!dst) {
        return;
    }
    memset(dst, 0, sizeof(*dst));
    dst->max_length = (size_t)-1;
    if (src) {
        *dst = *src;
        if (dst->max_length == 0) {
            dst->max_length = (size_t)-1;
        }
    }
}

static bool fc_iter_filter_dust_bam(
    bam1_t *b, const fc_iter_opts *opts, fc_file_stats *stats, dust_buffer *dust)
{
    if (!opts || !opts->dust) {
        return true;
    }

    size_t len = (size_t)b->core.l_qseq;
    uint8_t *seq_data = bam_get_seq(b);

    // Reuse/grow buffer as needed
    if (len + 1 > dust->size) {
        dust->buf = fc_realloc(dust->buf, len + 1, "bam dust seq");
        dust->size = len + 1;
    }
    char *seq = dust->buf;

    for (size_t i = 0; i < len; ++i) {
        seq[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    }
    seq[len] = '\0';

    bool pass = fc_filter_dust((uint8_t *)seq, len, opts->dust_t, opts->dust_w, opts->max_dust);
    if (!pass) {
        stats->failures[R_DUST_MASKED]++;
        return false;
    }
    return true;
}

static bool fc_iter_process_bam_record(
    fc_read *out, const bam_hdr_t *hdr, const kstring_t *orig_comment,
    const fc_iter_opts *opts, bool apply_filters, fc_file_stats *stats,
    double *q_sum, double *q_comp, dust_buffer *dust)
{
    size_t len = (size_t)out->bam->core.l_qseq;
    float mean_q = fc_bam_ensure_mean_quality(out->bam, opts->force_recalc_qual);

    if (apply_filters) {
        if (!fc_iter_filter_length(len, opts->min_length, opts->max_length, stats)) {
            return false;
        }
        if (!fc_iter_filter_quality(mean_q, opts->min_qscore, stats)) {
            return false;
        }
        if (!fc_iter_filter_dust_bam(out->bam, opts, stats, dust)) {
            return false;
        }
        if (opts->tag_name[0] != '\0') {
            if (!fc_filter_int_tag(out->bam, opts->tag_name, opts->tag_value, false)) {
                return false;
            }
        }
        if (opts->read_group != NULL) {
            if (!fc_filter_rg(out->bam, opts->read_group, false)) {
                return false;
            }
        }
    }

    fc_iter_accum_add(len, mean_q, stats, q_sum, q_comp);
    fc_read_reset_from_bam(out, hdr, orig_comment);
    return true;
}

static bool fastq_iter_process_record(fc_fastq_iter *iter, fc_read *out)
{
    return fc_iter_process_bam_record(out, NULL, &iter->seq->comment,
        &iter->opts, iter->apply_filters, &iter->stats,
        &iter->q_sum, &iter->q_comp, &iter->dust);
}

static bool bam_iter_process_record(fc_bam_iter *iter, fc_read *out)
{
    return fc_iter_process_bam_record(out, iter->hdr, NULL,
        &iter->opts, iter->apply_filters, &iter->stats,
        &iter->q_sum, &iter->q_comp, &iter->dust);
}

static int fc_iter_finish_common(
    fc_file_stats *stats,
    double q_sum,
    int status,
    bool add_status_failure,
    fc_file_stats *out_stats)
{
    if (!stats) {
        return EINVAL;
    }

    if (add_status_failure) {
        if (status < -1) {
            stats->failures[F_STREAM_ERROR]++;
        } else if (status == -1) {
            stats->failures[F_FILE_OK]++;
        } else if (status < 0) {
            stats->failures[F_UNKNOWN_ERROR]++;
        }
    }

    finalize_file_stats(stats, q_sum);

    if (out_stats) {
        *out_stats = *stats;
    }

    return (status == -1 || status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

static void fastq_iter_close_resources(fc_fastq_iter *iter)
{
    if (!iter) {
        return;
    }
    if (iter->seq) {
        kseq_destroy(iter->seq);
        iter->seq = NULL;
    }
    if (iter->fp) {
        gzclose(iter->fp);
        iter->fp = NULL;
    }
    free_dust_buffer(&iter->dust);
}

static bool iter_opts_enable_filters(const fc_iter_opts *opts)
{
    if (!opts) {
        return false;
    }
    return (
        opts->min_length > 0 || opts->max_length != (size_t)-1 || opts->min_qscore > 0.0f || opts->dust || opts->tag_name[0] != '\0' || opts->read_group != NULL);
}

static int fc_bam_iter_set_region(fc_bam_iter *iter, const char *chr, hts_pos_t start, hts_pos_t end);

static int fc_fastq_iter_open(fc_fastq_iter *iter, const char *fname, const fc_iter_opts *opts)
{
    if (!iter || !fname || fname[0] == '\0') {
        return EINVAL;
    }
    memset(iter, 0, sizeof(*iter));
    iter->fname = fname;
    fc_file_stats_init(&iter->stats);
    fc_iter_opts_init(&iter->opts, opts);
    iter->apply_filters = iter_opts_enable_filters(&iter->opts);

    if (strcmp(fname, "-") == 0) {
        int fd = dup(fileno(stdin));
        if (fd < 0) {
            FC_LOG_ERROR("could not duplicate stdin: %s", strerror(errno));
            iter->stats.failures[F_UNKNOWN_ERROR] = 1;
            return errno;
        }
        iter->fp = gzdopen(fd, "r");
        if (!iter->fp) {
            close(fd);
        }
    } else {
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
    }
    if (iter->fp == NULL) {
        FC_LOG_ERROR("could not open file %s: %s", fname, strerror(errno));
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        return errno;
    }

    iter->seq = kseq_init(iter->fp);
    if (!iter->seq) {
        FC_LOG_ERROR("could not allocate FASTQ iterator sequence buffer");
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        fastq_iter_close_resources(iter);
        return ENOMEM;
    }
    return 0;
}

static fc_iter_next_status fc_fastq_iter_next(fc_fastq_iter *iter, fc_read *out)
{
    if (!iter || !iter->seq || !out || !out->bam) {
        return FC_ITER_ERROR;
    }

    while ((iter->status = kseq_read(iter->seq)) != -1) {
        if (iter->status == -2) {
            iter->stats.failures[F_QUAL_TRUNCATED]++;
            iter->truncated = true;
            continue;
        } else if (iter->status < 0) { // -3 or other error
            iter->status = -3;
            return FC_ITER_ERROR;
        }
        if (iter->seq->qual.l == 0) {
            iter->stats.failures[F_QUAL_MISSING]++;
            iter->truncated = true;
            iter->status = -2;
            continue;
        }

        iter->truncated = false;
        iter->stats.failures[R_RECORD_OK]++;

        // convert from kseq to BAM record in-place for filtering and metadata extraction
        if (kseq_to_bam(out->bam, iter->seq) < 0) {
            FC_LOG_FATAL("Error building BAM record for read '%s'", iter->seq->name.s);
        }
        if (fastq_iter_process_record(iter, out)) {
            return FC_ITER_RECORD;
        }
    }

    if (iter->truncated) {
        iter->status = -3;
    } else {
        iter->status = -1;
    }

    return FC_ITER_EOF;
}

static int fc_fastq_iter_finish(fc_fastq_iter *iter, fc_file_stats *out_stats)
{
    if (!iter) {
        return EINVAL;
    }

    int rc = fc_iter_finish_common(&iter->stats, iter->q_sum, iter->status, true, out_stats);

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

    return rc;
}

static int fc_bam_iter_open(
    fc_bam_iter *iter,
    htsFile *fp,
    hts_idx_t *idx,
    sam_hdr_t *hdr,
    const fc_iter_opts *opts)
{
    if (!iter || !fp || !hdr) {
        return -1;
    }
    memset(iter, 0, sizeof(*iter));
    fc_file_stats_init(&iter->stats);
    fc_iter_opts_init(&iter->opts, opts);

    iter->fp = fp;
    iter->hdr = hdr;
    iter->idx = idx;
    iter->own_idx = false;

    const char *bed_path = iter->opts.bed;
    const char *region_string = iter->opts.region;
    iter->use_regions = (bed_path != NULL) || (region_string != NULL);
    iter->apply_filters = iter_opts_enable_filters(&iter->opts);
    if (iter->use_regions && iter->idx == NULL) {
        // When streaming regions without an index we must apply overlap filtering in fc_bam_iter_next().
        iter->apply_filters = true;
    }

    if (iter->use_regions) {
        if (iter->idx == NULL && fp->fn != NULL && strcmp(fp->fn, "-") != 0) {
            iter->idx = sam_index_load(fp, fp->fn);
            iter->own_idx = (iter->idx != NULL);
        }
        if (iter->idx == NULL && fp->fn != NULL && strcmp(fp->fn, "-") != 0) {
            FC_LOG_ERROR("Cannot find index for BAM when regions are requested.");
            return -1;
        }

        iter->reg = init_region_iterator(bed_path, region_string, hdr);

        regiter_status check = next_region(&iter->reg);
        if (check == REGITER_EOF) {
            iter->reg_active = false;
            return 0;
        }
        if (check == REGITER_PARSE_ERROR && bed_path == NULL && region_string != NULL) {
            FC_LOG_ERROR("Failed to parse region string: %s", region_string);
            return -1;
        }
        while (check != REGITER_OK) {
            check = next_region(&iter->reg);
            if (check == REGITER_EOF) {
                break;
            }
        }
        iter->reg_active = (check == REGITER_OK);
        if (iter->reg_active && iter->idx != NULL) {
            if (fc_bam_iter_set_region(iter, iter->reg.chr, iter->reg.start, iter->reg.end) != 0) {
                iter->reg_active = false;
            }
        }
        if (iter->idx == NULL) {
            iter->bam = fc_bam_iter_data_create(fp, idx, hdr, NULL, 0, INT64_MAX);
            if (iter->bam == NULL) {
                return -1;
            }
        }
        return 0;
    }

    iter->bam = fc_bam_iter_data_create(fp, idx, hdr, NULL, 0, INT64_MAX);
    if (iter->bam == NULL) {
        return -1;
    }

    return 0;
}

static int fc_bam_iter_set_region(fc_bam_iter *iter, const char *chr, hts_pos_t start, hts_pos_t end)
{
    if (iter->bam) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
    }
    iter->bam = fc_bam_iter_data_create(iter->fp, iter->idx, iter->hdr, chr, start, end);
    return iter->bam ? 0 : -1;
}

static fc_iter_next_status fc_bam_iter_next(fc_bam_iter *iter, fc_read *out)
{
    if (!iter || !out || !out->bam) {
        return FC_ITER_ERROR;
    }

    for (;;) {
        // If no active bam iterator, try to advance to next region (if any)
        if (!iter->bam) {
            if (!iter->use_regions || !iter->reg_active) {
                iter->status = -1;
                return FC_ITER_EOF;
            }
        }

        iter->status = iter->bam ? fc_bam_read(iter->bam, out->bam) : -1;

        if (iter->status < 0) {
            if (iter->use_regions) {
                // advance to next region
                regiter_status check = next_region(&iter->reg);
                while (check == REGITER_PARSE_ERROR || check == REGITER_REF_MISSING) { // skip invalid/unknown regions
                    check = next_region(&iter->reg);
                }
                if (check == REGITER_OK) {
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
                    iter->status = fc_bam_read(iter->bam, out->bam);
                    if (iter->status < 0) {
                        continue; // empty region, move on
                    }
                } else {
                    iter->reg_active = false;
                    iter->status = -1;
                    return FC_ITER_EOF;
                }
            } else {
                return (iter->status < -1) ? FC_ITER_ERROR : FC_ITER_EOF;
            }
        }

        if (iter->status < 0) {
            return (iter->status < -1) ? FC_ITER_ERROR : FC_ITER_EOF;
        }

        iter->stats.failures[R_RECORD_OK]++;

        if (iter->use_regions && iter->idx == NULL) {
            int tid = out->bam->core.tid;
            hts_pos_t pos = out->bam->core.pos;
            hts_pos_t rend = bam_endpos(out->bam);
            while (iter->reg_active && (tid > iter->reg.tid || (tid == iter->reg.tid && pos >= iter->reg.end))) {
                regiter_status check = next_region(&iter->reg);
                if (check != REGITER_OK) {
                    iter->reg_active = false;
                    break;
                }
            }
            if (!iter->reg_active) {
                iter->stats.failures[R_RECORD_OK]--;
                iter->status = -1;
                return FC_ITER_EOF;
            }
            bool overlap = (tid == iter->reg.tid) && (pos < iter->reg.end) && (rend > iter->reg.start);
            if (!overlap) {
                iter->stats.failures[R_RECORD_OK]--;
                continue;
            }
        }

        if (bam_iter_process_record(iter, out)) {
            return FC_ITER_RECORD;
        }
    }

    return FC_ITER_EOF;
}

static int fc_bam_iter_finish(fc_bam_iter *iter, fc_file_stats *out_stats)
{
    if (!iter) {
        return EINVAL;
    }
    return fc_iter_finish_common(&iter->stats, iter->q_sum, iter->status, true, out_stats);
}

static void fc_bam_iter_close(fc_bam_iter *iter)
{
    if (!iter) {
        return;
    }
    if (iter->bam) {
        fc_bam_iter_data_destroy(iter->bam);
        iter->bam = NULL;
    }
    if (iter->hdr) {
        sam_hdr_destroy(iter->hdr);
        iter->hdr = NULL;
    }
    if (iter->fp) {
        hts_close(iter->fp);
        iter->fp = NULL;
    }
    if (iter->use_regions) {
        destroy_region_iterator(&iter->reg);
        memset(&iter->reg, 0, sizeof(iter->reg));
        iter->use_regions = false;
    }
    if (iter->own_idx && iter->idx) {
        hts_idx_destroy(iter->idx);
        iter->idx = NULL;
    }
    iter->own_idx = false;
    free_dust_buffer(&iter->dust);
}

static int fc_bam_iter_open_path(fc_bam_iter *bmi, const char *path, const fc_iter_opts *opts)
{
    if (!bmi || !path || path[0] == '\0') {
        return EINVAL;
    }
    fc_file_stats_init(&bmi->stats);

    htsFile *bam_fp = hts_open(path, "rb");
    if (!bam_fp && strcmp(path, "-") == 0) {
        bam_fp = hts_open(path, "r");
    }
    if (!bam_fp) {
        bmi->stats.failures[F_UNKNOWN_ERROR] = 1;
        return errno ? errno : EIO;
    }

    // attach caller-provided thread pool, if any
    if (opts && opts->thread_pool) {
        hts_set_opt(bam_fp, HTS_OPT_THREAD_POOL, opts->thread_pool);
    }

    sam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);
    if (!bam_hdr) {
        hts_close(bam_fp);
        bmi->stats.failures[F_UNKNOWN_ERROR] = 1;
        return EINVAL;
    }

    int rc = fc_bam_iter_open(bmi, bam_fp, NULL, bam_hdr, opts);
    if (rc != 0) {
        fc_bam_iter_close(bmi);
        bmi->stats.failures[F_UNKNOWN_ERROR] = 1;
        return rc;
    }
    return 0;
}

static bool has_suffix(const char *path, const char *suffix)
{
    size_t lp = strlen(path);
    size_t ls = strlen(suffix);
    if (lp < ls) {
        return false;
    }
    return strcmp(path + (lp - ls), suffix) == 0;
}

bool fc_path_matches_kind(const char *path, fc_file_kind expected_kind)
{
    if (!path) {
        return false;
    }

    if (expected_kind == FC_FILE_FASTQ) {
        return has_suffix(path, ".fastq") || has_suffix(path, ".fastq.gz") || has_suffix(path, ".fq") || has_suffix(path, ".fq.gz");
    } else if (expected_kind == FC_FILE_BAM) {
        return has_suffix(path, ".bam") || has_suffix(path, ".sam");
    }
    return false;
}

static int append_file(char ***files, size_t *n_files, const char *path)
{
    size_t new_count = *n_files + 1;
    *files = fc_realloc(*files, new_count * sizeof(**files), "file list array");

    char *dup = fc_strdup(path, "path string for region parsing");
    (*files)[*n_files] = dup;
    *n_files = new_count;
    return 0;
}

typedef struct {
    char *path;
    int depth;
} pending_path;

static int collect_files_from_roots(
    const fc_path_list *input_paths, int recurse,
    fc_file_kind expected_kind, fc_path_list *out_paths)
{
    int status = 0;
    char **files = NULL;
    size_t n_files = 0;

    pending_path *queue = NULL;
    size_t qlen = 0, qcap = 0, qhead = 0;

    if (!input_paths || !input_paths->paths || input_paths->n_paths == 0) {
        return EINVAL;
    }

    // seed queue
    for (size_t i = 0; i < input_paths->n_paths; ++i) {
        if (qlen == qcap) {
            size_t newcap = qcap ? qcap * 2 : 8;
            queue = fc_realloc(queue, newcap * sizeof(*queue), "path iteration queue");
            qcap = newcap;
        }
        queue[qlen].path = fc_strdup(input_paths->paths[i], "queued path for iteration");
        queue[qlen].depth = recurse;
        qlen++;
    }

    while (qhead < qlen) {
        pending_path cur = queue[qhead++];
        if (strcmp(cur.path, "-") == 0) {
            int rc = append_file(&files, &n_files, cur.path);
            if (rc != 0) {
                status = (status == 0) ? rc : status;
                free(cur.path);
                break;
            }
            free(cur.path);
            continue;
        }
        struct stat st = { 0 };
        if (stat(cur.path, &st) != 0) {
            FC_LOG_ERROR("could not stat '%s': %s", cur.path, strerror(errno));
            status = (status == 0) ? errno : status;
            free(cur.path);
            continue;
        }
        if (S_ISDIR(st.st_mode)) {
            // depth semantics:
            //   0  => do not descend
            //   N>0 => descend up to N levels (decrement each level)
            //   N<0 => unbounded recursion
            if (cur.depth == 0) {
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
                if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
                    continue;
                }
                char *child = fc_asprintf("recursive input child path", "%s/%s", cur.path, ent->d_name);
                if (qlen == qcap) {
                    size_t newcap = qcap ? qcap * 2 : 8;
                    queue = fc_realloc(queue, newcap * sizeof(*queue), "path iteration queue (recursive)");
                    qcap = newcap;
                }
                queue[qlen].path = child;
                queue[qlen].depth = (cur.depth < 0) ? -1 : (cur.depth - 1);
                qlen++;
            }
            closedir(dir);
            free(cur.path);
            continue;
        }

        bool matches_expected = fc_path_matches_kind(cur.path, expected_kind);
        bool matches_other = fc_path_matches_kind(
            cur.path, expected_kind == FC_FILE_FASTQ ? FC_FILE_BAM : FC_FILE_FASTQ);

        if (matches_expected) {
            int rc = append_file(&files, &n_files, cur.path);
            if (rc != 0) {
                status = (status == 0) ? rc : status;
                free(cur.path);
                break;
            }
            free(cur.path);
            continue;
        }
        if (matches_other) {
            FC_LOG_ERROR("Input '%s' does not match required input type for this command.", cur.path);
            status = (status == 0) ? EINVAL : status;
        }
        // Unknown extension: ignore
        free(cur.path);
    }

    for (size_t i = qhead; i < qlen; ++i) {
        free(queue[i].path);
    }
    free(queue);
    if (status != 0) {
        for (size_t i = 0; i < n_files; ++i) {
            free(files[i]);
        }
        free(files);
        files = NULL;
        n_files = 0;
    }
    out_paths->paths = files;
    out_paths->n_paths = n_files;
    return status;
}

int fc_iter_collect_paths(
    const fc_path_list *input_paths, int recurse,
    fc_file_kind expected_kind, bool stdin_is_path_list, fc_path_list *out_paths)
{
    if (!out_paths) {
        return EINVAL;
    }

    int rc = 0;
    out_paths->paths = NULL;
    out_paths->n_paths = 0;

    if (stdin_is_path_list) {
        rc = fc_read_paths_from_fofn_sources(input_paths, expected_kind, out_paths);
    } else if (input_paths && input_paths->n_paths == 1 && input_paths->paths && strcmp(input_paths->paths[0], "-") == 0) {
        // stdin carries data stream (FASTQ/BAM), not a path list
        char **files = fc_alloc(1, sizeof(*files), "stdin path array");
        files[0] = fc_strdup("-", "stdin path marker");
        out_paths->paths = files;
        out_paths->n_paths = 1;
    } else {
        rc = collect_files_from_roots(input_paths, recurse, expected_kind, out_paths);
    }

    return rc;
}

int fc_iter_init(
    fc_iter *it, const fc_path_list *input_paths,
    int recurse, bool verbose, fc_file_kind expected_kind, bool stdin_is_path_list)
{
    if (!it || !input_paths || !input_paths->paths || input_paths->n_paths == 0) {
        return EINVAL;
    }

    int rc = 0;
    fc_path_list collected = { 0 };
    memset(it, 0, sizeof(*it));
    it->recurse = recurse;
    it->verbose = verbose;
    it->expected_kind = expected_kind;
    rc = fc_iter_collect_paths(
        input_paths, recurse, expected_kind, stdin_is_path_list, &collected);
    if (rc == 0) {
        it->files = collected.paths;
        it->n_files = collected.n_paths;
    }

    return rc;
}

bool fc_iter_next_file(fc_iter *it, const char **out_file)
{
    if (!it || it->idx >= it->n_files) {
        return false;
    }

    it->current = it->files[it->idx++];
    if (it->verbose) {
        FC_LOG_INFO("Processing %s", it->current);
    }
    if (out_file) {
        *out_file = it->current;
    }

    return true;
}

int fc_iter_open_record_cursor(
    fc_iter *it,
    const fc_iter_opts *opts)
{
    if (!it || !it->current) {
        return EINVAL;
    }

    int rc = EINVAL;
    it->cursor_open = false;
    if (it->expected_kind == FC_FILE_FASTQ) {
        // When iterating multiple FASTQ inputs, release any previous FASTQ cursor
        // before opening the next file.
        fastq_iter_close_resources(&it->fqi);
        fc_file_stats_init(&it->fqi.stats);
        rc = fc_fastq_iter_open(&it->fqi, it->current, opts);
    } else {
        // BAM path
        // When iterating multiple BAM inputs (fastcat fastq --input-fmt bam),
        // release any previous BAM cursor before opening the next file.
        fc_bam_iter_close(&it->bmi);
        memset(&it->bmi, 0, sizeof(it->bmi));
        fc_file_stats_init(&it->bmi.stats);
        rc = fc_bam_iter_open_path(&it->bmi, it->current, opts);
    }
    it->cursor_open = (rc == 0);

    return rc;
}

fc_iter_next_status fc_iter_next_record(fc_iter *it, fc_read *out)
{
    fc_iter_next_status status = FC_ITER_ERROR;
    if (!it || !it->current || !out) {
        return status;
    }

    if (it->expected_kind == FC_FILE_FASTQ) {
        status = fc_fastq_iter_next(&it->fqi, out);
    } else {
        status = fc_bam_iter_next(&it->bmi, out);
    }
    return status;
}

int fc_iter_finish_file(fc_iter *it)
{
    if (!it || !it->current) {
        return EINVAL;
    }

    int rc = EINVAL;
    if (it->expected_kind == FC_FILE_FASTQ) {
        rc = fc_fastq_iter_finish(&it->fqi, NULL);
    } else {
        rc = fc_bam_iter_finish(&it->bmi, NULL);
    }
    it->cursor_open = false;

    return rc;
}

void fc_iter_close(fc_iter *it)
{
    if (!it) {
        return;
    }

    if (it->expected_kind == FC_FILE_FASTQ) {
        fastq_iter_close_resources(&it->fqi);
    } else {
        fc_bam_iter_close(&it->bmi);
    }

    fc_path_list paths = { .paths = it->files, .n_paths = it->n_files };
    fc_path_list_free(&paths);
    it->files = paths.paths;
    it->n_files = paths.n_paths;
    it->idx = 0;
    it->current = NULL;
}

void fc_iter_merge_file_stats(const fc_iter *it, fc_file_stats *dst)
{
    if (!it || !dst) {
        return;
    }

    const fc_file_stats *src = fc_iter_stats(it);
    if (!src) {
        return;
    }

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

const fc_file_stats *fc_iter_stats(const fc_iter *it)
{
    if (!it || !it->current) {
        return NULL;
    }

    return it->expected_kind == FC_FILE_FASTQ ? &it->fqi.stats : &it->bmi.stats;
}

void fc_path_list_free(fc_path_list *paths)
{
    if (!paths || !paths->paths) {
        return;
    }

    for (size_t i = 0; i < paths->n_paths; ++i) {
        free(paths->paths[i]);
    }
    free(paths->paths);
    paths->paths = NULL;
    paths->n_paths = 0;
}
