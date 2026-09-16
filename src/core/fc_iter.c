
#include "core/fc_iter.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "htslib/bgzf.h"

#include "core/fc_bam_utils.h"
#include "core/fc_filter.h"
#include "core/fc_log.h"
#include "core/fc_run_summary.h"
#include "core/kseq_gz.h"

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

void fc_iter_opts_disable_filters(fc_iter_opts *opts)
{
    if (opts == NULL) {
        return;
    }
    opts->min_length = 0;
    opts->max_length = (size_t)-1;
    opts->min_qscore = 0.0f;
    opts->dust = false;
    opts->read_group = NULL;
    opts->tag_name[0] = '\0';
    opts->tag_name[1] = '\0';
    opts->tag_value = -1;
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
    fc_bam_rg_cache *rg_cache,
    const fc_iter_opts *opts, bool apply_filters, fc_file_stats *stats,
    double *q_sum, double *q_comp, dust_buffer *dust, size_t barcode)
{
    size_t len = (size_t)out->bam->core.l_qseq;
    float mean_q = fc_bam_ensure_mean_quality(out->bam, opts->force_recalc_qual);
    // non-finite qs/QS values are treated as low quality
    // we've already filtered out reads with missing quality scores
    if (!isfinite(mean_q)) {
        stats->failures[R_LOW_QUALITY]++;
        fc_run_summary_add_filter_event(
            opts->run_summary, opts->use_barcodes ? barcode : 0, R_LOW_QUALITY);
        return false;
    }

    if (apply_filters) {
        if (!fc_iter_filter_length(len, opts->min_length, opts->max_length, stats)) {
            fc_run_summary_add_filter_event(opts->run_summary, opts->use_barcodes ? barcode : 0,
                len > opts->max_length ? R_TOO_LONG : R_TOO_SHORT);
            return false;
        }
        if (!fc_iter_filter_quality(mean_q, opts->min_qscore, stats)) {
            fc_run_summary_add_filter_event(
                opts->run_summary, opts->use_barcodes ? barcode : 0, R_LOW_QUALITY);
            return false;
        }
        if (!fc_iter_filter_dust_bam(out->bam, opts, stats, dust)) {
            fc_run_summary_add_filter_event(
                opts->run_summary, opts->use_barcodes ? barcode : 0, R_DUST_MASKED);
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
    fc_read_reset_from_bam(out, hdr, rg_cache, orig_comment);
    return true;
}

static bool fastq_iter_process_record(fc_fastq_iter *iter, fc_read *out,
    const kstring_t *comment, const fc_fastq_record *parsed_record)
{
    size_t barcode = fc_read_barcode_index_from_bam(out->bam, NULL, NULL);
    fc_run_summary_add_filter_event(
        iter->opts.run_summary, iter->opts.use_barcodes ? barcode : 0, R_RECORD_SEEN);
    bool accepted = fc_iter_process_bam_record(out, NULL, comment, iter->rg_cache,
        &iter->opts, iter->apply_filters, &iter->stats,
        &iter->q_sum, &iter->q_comp, &iter->dust, barcode);
    if (accepted) {
        // kseq storage remains valid until the next fc_iter_next_record call.
        // The writer consumes this record before then, allowing FASTQ output
        // to avoid unpacking SEQ and QUAL from its BAM representation.
        out->fastq = parsed_record ? NULL : iter->seq;
        out->fastq_record = parsed_record;
        fc_run_summary_add_filter_event(
            iter->opts.run_summary, iter->opts.use_barcodes ? barcode : 0, R_RECORD_OK);
    }
    return accepted;
}

static bool bam_iter_process_record(fc_bam_iter *iter, fc_read *out)
{
    size_t barcode = fc_read_barcode_index_from_bam(out->bam, iter->hdr, iter->rg_cache);
    fc_run_summary_add_filter_event(
        iter->opts.run_summary, iter->opts.use_barcodes ? barcode : 0, R_RECORD_SEEN);
    bool accepted = fc_iter_process_bam_record(out, iter->hdr, NULL, iter->rg_cache,
        &iter->opts, iter->apply_filters, &iter->stats,
        &iter->q_sum, &iter->q_comp, &iter->dust, barcode);
    if (accepted) {
        fc_run_summary_add_filter_event(
            iter->opts.run_summary, iter->opts.use_barcodes ? barcode : 0, R_RECORD_OK);
    }
    return accepted;
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
    if (iter->parser) {
        fc_fastq_parser_destroy(iter->parser);
        iter->parser = NULL;
    }
    if (iter->seq) {
        kseq_destroy(iter->seq);
        iter->seq = NULL;
    }
    if (iter->fp) {
        (void)fc_fastq_stream_close(iter->fp);
        iter->fp = NULL;
    }
    if (iter->rg_cache) {
        fc_bam_rg_cache_destroy(iter->rg_cache);
        iter->rg_cache = NULL;
    }
    free_dust_buffer(&iter->dust);
}

static bool bam_iter_input_offset(const fc_bam_iter *iter, uint64_t *out)
{
    if (!iter || !iter->fp || !out || iter->fp->format.format != bam
        || (iter->use_regions && iter->idx != NULL)) {
        return false;
    }
    BGZF *bgzf = iter->fp->fp.bgzf;
    if (!bgzf) {
        return false;
    }
    int64_t virtual_offset = bgzf_tell(bgzf);
    if (virtual_offset < 0) {
        return false;
    }
    *out = (uint64_t)virtual_offset >> 16;
    return true;
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
        iter->fp = fc_fastq_stream_open_fd(fd, iter->opts.fastq_threads);
        if (!iter->fp) {
            close(fd);
        }
    } else {
        fc_file_stat finfo;
        if (fc_stat_path(fname, &finfo) == -1) {
            FC_LOG_ERROR("could not process file %s: %s", fname, strerror(errno));
            iter->stats.failures[F_UNKNOWN_ERROR] = 1;
            return errno;
        }

        if ((finfo.st_mode & S_IFMT) == S_IFDIR) {
            iter->stats.failures[F_UNKNOWN_ERROR] = 1;
            return EISDIR;
        }

        iter->fp = fc_fastq_stream_open_path(fname, iter->opts.fastq_threads);
    }
    if (iter->fp == NULL) {
        FC_LOG_ERROR("could not open file %s: %s", fname, strerror(errno));
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        return errno;
    }

#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    iter->parser = fc_fastq_parser_create(iter->fp);
    if (!iter->parser) {
        FC_LOG_ERROR("could not allocate FASTQ parser");
#else
    iter->seq = kseq_init(iter->fp);
    if (!iter->seq) {
        FC_LOG_ERROR("could not allocate FASTQ iterator sequence buffer");
#endif
        iter->stats.failures[F_UNKNOWN_ERROR] = 1;
        fastq_iter_close_resources(iter);
        return ENOMEM;
    }
    iter->rg_cache = fc_bam_rg_cache_create(NULL);
    if (NULL == iter->rg_cache) {
        fastq_iter_close_resources(iter);
        return ENOMEM;
    }
    return 0;
}

static fc_iter_next_status fc_fastq_iter_next(fc_fastq_iter *iter, fc_read *out)
{
    if (!iter || (!iter->parser && !iter->seq) || !out || !out->bam) {
        return FC_ITER_ERROR;
    }

#if !defined(NOTHREADS) && !defined(NOPURRFLATE)
    for (;;) {
        int status = fc_fastq_parser_next(iter->parser, out->bam, &iter->fastq_record);
        if (status == 0) {
            iter->status = iter->truncated ? -3 : -1;
            return FC_ITER_EOF;
        }
        if (status == -2) {
            iter->stats.failures[F_QUAL_TRUNCATED]++;
            iter->truncated = true;
            continue;
        }
        if (status == -4) {
            iter->stats.failures[F_QUAL_MISSING]++;
            iter->truncated = true;
            continue;
        }
        if (status < 0) {
            iter->status = -3;
            return FC_ITER_ERROR;
        }

        iter->truncated = false;
        iter->records_read++;
        iter->bases_read += iter->fastq_record.sequence_len;
        iter->stats.failures[R_RECORD_SEEN]++;

        kstring_t comment = {
            .l = iter->fastq_record.comment_len,
            .m = iter->fastq_record.comment_len,
            .s = (char *)iter->fastq_record.comment,
        };
        if (fastq_iter_process_record(
                iter, out, &comment, &iter->fastq_record)) {
            iter->stats.failures[R_RECORD_OK]++;
            return FC_ITER_RECORD;
        }
    }
#else
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
        iter->records_read++;
        iter->bases_read += iter->seq->seq.l;
        iter->stats.failures[R_RECORD_SEEN]++;

        if (fc_fastq_record_to_bam(out->bam,
                iter->seq->name.s, iter->seq->name.l,
                iter->seq->comment.s, iter->seq->comment.l,
                iter->seq->seq.s, iter->seq->qual.s, iter->seq->seq.l)
            < 0) {
            FC_LOG_FATAL("Error building BAM record for read '%s'", iter->seq->name.s);
        }
        if (fastq_iter_process_record(iter, out, &iter->seq->comment, NULL)) {
            iter->stats.failures[R_RECORD_OK]++;
            return FC_ITER_RECORD;
        }
    }

    if (iter->truncated) {
        iter->status = -3;
    } else {
        iter->status = -1;
    }

    return FC_ITER_EOF;
#endif
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
    iter->rg_cache = fc_bam_rg_cache_create(hdr);
    if (NULL == iter->rg_cache) {
        return -1;
    }
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
        if (iter->reg.error != 0) {
            return -1;
        }

        regiter_status check = next_region(&iter->reg);
        if (check == REGITER_EOF) {
            iter->reg_active = false;
            return 0;
        }
        // Tolerate malformed BED entries and continue with subsequent lines.
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
                // A stream is scanned once and filtered against the complete
                // normalised region set below; it cannot seek to another region.
                if (iter->idx == NULL) {
                    return (iter->status < -1) ? FC_ITER_ERROR : FC_ITER_EOF;
                }
                // advance to next region
                regiter_status check = next_region(&iter->reg);
                // Tolerate malformed entries and unknown references instead
                // of aborting the run.
                while (check == REGITER_PARSE_ERROR || check == REGITER_REF_MISSING) {
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

        if (iter->use_regions && iter->idx == NULL) {
            int tid = out->bam->core.tid;
            hts_pos_t pos = out->bam->core.pos;
            hts_pos_t rend = bam_endpos(out->bam);
            if (!fc_bam_coordinate_order_check(&iter->coordinate_order, out->bam)) {
                FC_LOG_ERROR("Region selection from an unindexed BAM stream requires coordinate-sorted records.");
                iter->hard_error = true;
                return FC_ITER_ERROR;
            }
            while (iter->reg_active && (tid > iter->reg.tid || (tid == iter->reg.tid && pos >= iter->reg.end))) {
                regiter_status check = next_region(&iter->reg);
                // Skip malformed entries like indexed mode does.
                while (check == REGITER_PARSE_ERROR || check == REGITER_REF_MISSING) {
                    check = next_region(&iter->reg);
                }
                if (check != REGITER_OK) {
                    iter->reg_active = false;
                    break;
                }
            }
            if (!iter->reg_active) {
                iter->status = -1;
                return FC_ITER_EOF;
            }
            bool overlap = (tid == iter->reg.tid) && (pos < iter->reg.end) && (rend > iter->reg.start);
            if (!overlap) {
                continue;
            }
        }

        iter->records_read++;
        iter->bases_read += out->bam->core.l_qseq;

        if (fc_bam_qualities_missing(out->bam)) {
            if (fc_is_canonical_read_record(out->bam)) {
                iter->stats.failures[F_QUAL_MISSING]++;
            }
            // QUAL=* is valid for records such as secondary alignments. Let
            // the writer record its alignment flags. Missing quality on a
            // canonical read is tracked as a failure; additional alignments
            // are valid records but never enter read-level models or sinks.
            fc_read_reset_from_bam(out, iter->hdr, iter->rg_cache, NULL);
            return FC_ITER_RECORD;
        }

        iter->stats.failures[R_RECORD_SEEN]++;
        if (bam_iter_process_record(iter, out)) {
            iter->stats.failures[R_RECORD_OK]++;
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
    if (iter->rg_cache) {
        fc_bam_rg_cache_destroy(iter->rg_cache);
        iter->rg_cache = NULL;
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

int fc_iter_init(
    fc_iter *it, fc_resolved_inputs *inputs, bool verbose)
{
    if (!it || !inputs || !inputs->paths.paths || inputs->paths.n_paths == 0) {
        return EINVAL;
    }
    memset(it, 0, sizeof(*it));
    it->verbose = verbose;
    it->expected_kind = inputs->kind;
    it->files = inputs->paths.paths;
    it->n_files = inputs->paths.n_paths;
    inputs->paths.paths = NULL;
    inputs->paths.n_paths = 0;

    it->file_sizes = fc_alloc(it->n_files, sizeof(*it->file_sizes), "input file sizes");
    it->total_bytes_known = true;
    for (size_t i = 0; i < it->n_files; ++i) {
        fc_file_stat file_info;
        if (strcmp(it->files[i], "-") == 0
            || fc_stat_path(it->files[i], &file_info) != 0
            || !S_ISREG(file_info.st_mode)
            || file_info.st_size < 0
            || (uint64_t)file_info.st_size > UINT64_MAX - it->total_bytes) {
            it->total_bytes_known = false;
            continue;
        }
        it->file_sizes[i] = (uint64_t)file_info.st_size;
        it->total_bytes += it->file_sizes[i];
    }
    return 0;
}

bool fc_iter_next_file(fc_iter *it, const char **out_file)
{
    if (!it || it->idx >= it->n_files) {
        return false;
    }

    it->current_index = it->idx;
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
        // When iterating multiple BAM inputs, release any previous BAM cursor
        // before opening the next file.
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

    uint64_t current_bytes = 0;
    bool current_bytes_known = it->expected_kind == FC_FILE_FASTQ
        ? fc_fastq_stream_input_offset(it->fqi.fp, &current_bytes)
        : bam_iter_input_offset(&it->bmi, &current_bytes);
    uint64_t records_read = it->expected_kind == FC_FILE_FASTQ
        ? it->fqi.records_read
        : it->bmi.records_read;
    uint64_t bases_read = it->expected_kind == FC_FILE_FASTQ
        ? it->fqi.bases_read
        : it->bmi.bases_read;

    int rc = EINVAL;
    if (it->expected_kind == FC_FILE_FASTQ) {
        rc = fc_fastq_iter_finish(&it->fqi, NULL);
    } else {
        rc = fc_bam_iter_finish(&it->bmi, NULL);
    }
    it->cursor_open = false;
    it->completed_records += records_read;
    it->completed_bases += bases_read;
    if (!current_bytes_known) {
        it->total_bytes_known = false;
    } else if (it->total_bytes_known && rc == EXIT_SUCCESS) {
        it->completed_bytes += it->file_sizes[it->current_index];
    } else {
        it->completed_bytes += min(current_bytes, it->file_sizes[it->current_index]);
    }

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
    free(it->file_sizes);
    it->file_sizes = NULL;
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

bool fc_iter_has_hard_error(const fc_iter *it)
{
    return it != NULL && it->expected_kind == FC_FILE_BAM && it->bmi.hard_error;
}

void fc_iter_get_progress(const fc_iter *it, fc_iter_progress *out)
{
    if (!out) {
        return;
    }
    memset(out, 0, sizeof(*out));
    if (!it) {
        return;
    }

    out->bytes_read = it->completed_bytes;
    out->total_bytes = it->total_bytes;
    out->records_read = it->completed_records;
    out->bases_read = it->completed_bases;
    out->total_bytes_known = it->total_bytes_known;
    if (!it->cursor_open) {
        return;
    }

    uint64_t current_bytes = 0;
    bool current_bytes_known = it->expected_kind == FC_FILE_FASTQ
        ? fc_fastq_stream_input_offset(it->fqi.fp, &current_bytes)
        : bam_iter_input_offset(&it->bmi, &current_bytes);
    if (!current_bytes_known) {
        out->total_bytes_known = false;
    } else {
        out->bytes_read += min(current_bytes, it->file_sizes[it->current_index]);
    }
    if (it->expected_kind == FC_FILE_FASTQ) {
        out->records_read += it->fqi.records_read;
        out->bases_read += it->fqi.bases_read;
    } else {
        out->records_read += it->bmi.records_read;
        out->bases_read += it->bmi.bases_read;
    }
}
