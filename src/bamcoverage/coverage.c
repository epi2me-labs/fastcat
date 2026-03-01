#include "coverage.h"

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "../fc_log.h"
#include "common.h"
#include "regiter.h"

static inline int is_ref_and_query(const uint32_t op)
{
    // NOTE: CIGAR mode excludes deletions (BAM_CDEL) from coverage. This is
    // intentionally different from non-CIGAR mode, which uses POS..bam_endpos()
    // and therefore counts deleted/skipped reference positions as covered span.
    // This is the same as mosdepth's default behavior; our spanning mode is
    // mosdepth's --fast-mode.

    // RNA-seq workflows might want to count deletions while excluding introns
    // (BAM_CREF_SKIP), so we should maybe add a mode to CLI for configuring
    // whether either/or DEL and SKIP count as coverage.
    const int code = bam_cigar_op(op);
    return code == BAM_CMATCH || code == BAM_CEQUAL || code == BAM_CDIFF;
}

static inline int consumes_ref_only(const uint32_t op)
{
    const int code = bam_cigar_op(op);
    return code == BAM_CDEL || code == BAM_CREF_SKIP;
}

static inline void _write_bedgraph_line(
    BGZF *fp, const char *chrom, int64_t s, int64_t e, uint32_t cov)
{
    if (NULL == fp) {
        return;
    }
    //  (chrom + tabs + ints) easily fits in <128 bytes
    char buf[160];
    FC_SNPRINTF_CHECK(buf, sizeof(buf), "%s\t%" PRId64 "\t%" PRId64 "\t%u\n",
        chrom, s, e, cov);
    size_t n = strlen(buf);
    if (bgzf_write_small(fp, buf, n) != (ssize_t)n) {
        FC_LOG_FATAL("bgzf_write failed");
    }
}

static void coverage_fail_invalid_alignment(cov_writer w, const bam1_t *b, const char *reason)
{
    const char *chrom = "<unknown>";
    if (w && w->hdr && b && b->core.tid >= 0 && b->core.tid < w->hdr->n_targets) {
        chrom = w->hdr->target_name[b->core.tid];
    }
    FC_LOG_FATAL(
        "Invalid alignment coordinates for read '%s' on '%s': %s.",
        b ? bam_get_qname(b) : "<unknown>", chrom, reason);
}

static inline void _write_summary(
    FILE *fh_dist, FILE *fh_summary, const bed_region reg, const int64_t *stats, const int64_t *dist, size_t max_cover, uint32_t *thresholds, size_t n_thresholds)
{
    if (fh_dist == NULL && fh_summary == NULL) {
        return; // nothing to write
    }

    // write all distribution file entries and collect threshold values
    // for summary file.
    int64_t *thresh_counts = fc_alloc(n_thresholds, sizeof(int64_t), "threshold counts");
    ssize_t cur_thresh = n_thresholds - 1;
    double cum = 0.0;
    double positions = (double)stats[3];
    for (ssize_t cover = max_cover - 1; cover >= 0; --cover) {
        cum += (double)dist[cover];

        // dense CDF - mosdepth calls this "distribution"
        if (fh_dist) {
            double frac = cum / positions;
            if (frac < 1e-3) {
                continue; // mosdepth uses 8e-5 \:D/
            }
            fprintf(fh_dist, "%s\t%zu\t%.3f\n", reg->chr, cover, frac); // fh_dist only written for totals
        }

        // sparse CDF - mosdepth calls this "thresholds"
        // NOTE: we do this in long format not wide like mosdepth
        if (cover == thresholds[cur_thresh]) {
            thresh_counts[cur_thresh] = (int64_t)cum;
            if (cur_thresh > 0) {
                --cur_thresh;
            }
        }
    }

    // write summary file
    double mean = stats[3] == 0 ? 0 : (double)stats[2] / stats[3];
    int64_t min = stats[0] == INT64_MAX ? 0 : stats[0];
    fprintf(
        fh_summary, "%s\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%.2f\t%" PRId64 "\t%" PRId64,
        reg->chr, reg->start, reg->end, stats[3], stats[2], mean, min, stats[1]);
    for (size_t i = 0; i < n_thresholds; ++i) {
        double frac = stats[3] == 0 ? 0 : (double)thresh_counts[i] / stats[3];
        fprintf(fh_summary, "\t%.3f", frac);
    }
    fprintf(fh_summary, "\n");

    free(thresh_counts);
}

static void _fill_skipped_regions(cov_writer w)
{
    if (w == NULL || w->tid < 0) {
        return;
    }

    // We want to enforce explicit zeros for all regions in the input
    // BED files. The BED files are sorted (see regiter.c) so we can
    // detect anything we'd otherwise skip over.
    for (size_t i = 0; i < w->n_beds; ++i) {
        cov_writer_region wr = w->writers[i];
        for (size_t j = wr->cur_region; j < wr->bed->n_regions; ++j) {
            bed_region reg = &wr->bed->regions[j];
            char *region_name = wr->whole_chrom ? reg->chr : region_to_string(reg);
            if (reg->tid >= w->tid) {
                break;
            }
            int64_t stats[4] = { 0, 0, 0, reg->end - reg->start };
            int64_t *dist = fc_alloc(wr->max_cover, sizeof(int64_t), "coverage region distribution");
            _write_summary(wr->fh_dist, wr->fh_summary, reg, stats,
                dist, wr->max_cover, wr->thresholds, wr->n_thresholds);
            _write_bedgraph_line(wr->fh_fwd, reg->chr, reg->start, reg->end, 0);
            _write_bedgraph_line(wr->fh_rev, reg->chr, reg->start, reg->end, 0);
            _write_bedgraph_line(wr->fh_tot, reg->chr, reg->start, reg->end, 0);
            wr->cur_region++;
            // for completeness to avoid special casing, trivial update to stats here
            wr->stats[0] = 0;
            wr->stats[3] += stats[3];
            wr->dist[0] += (int64_t)(reg->end - reg->start);
            free(dist);
            if (!wr->whole_chrom) {
                free(region_name);
            }
        }
    }
}

static void _flush_contig(cov_writer w)
{
    if (w == NULL || w->tid < 0) {
        return;
    }

    // flush any regions in BED that won't otherwise get picked up
    _fill_skipped_regions(w);

    // first do a cumulative sum of deltas, this gives us the coverage trace
    for (int64_t pos = 1; pos < w->contig_len; ++pos) {
        w->diff_fwd[pos] += w->diff_fwd[pos - 1];
        w->diff_rev[pos] += w->diff_rev[pos - 1];
    }

    // for each BED file (or whole genome)
    for (size_t i = 0; i < w->n_beds; ++i) {
        cov_writer_region wr = w->writers[i];

        // for each region in BED file
        for (size_t j = wr->cur_region; j < wr->bed->n_regions; ++j) {
            bed_region reg = &wr->bed->regions[j];

            if (reg->tid != w->tid) {
                // moved to a new contig, stop processing
                wr->cur_region = j;
                break;
            }
            if (reg->start >= w->contig_len || reg->end <= 0) {
                wr->cur_region = j + 1;
                continue;
            }
            char *region_name = wr->whole_chrom ? reg->chr : region_to_string(reg);

            // chromosome summary
            int64_t stats[4] = { INT64_MAX, 0, 0, reg->end - reg->start }; // min, max, total, positions
            size_t max_cover = wr->max_cover; // which goes up to all thresholds
            int64_t *dist = fc_alloc(max_cover, sizeof(int64_t), "coverage contig distribution");
            for (int64_t pos = reg->start; pos < reg->end; ++pos) {
                int64_t cov = w->diff_fwd[pos] + w->diff_rev[pos];
                // Negative coverage should be impossible given the diff array logic.
                // This assert is primarily for sanitizer validation during development.
                assert(cov >= 0 && "Coverage should never be negative");
                stats[0] = min(stats[0], cov);
                stats[1] = max(stats[1], cov);
                stats[2] += cov;
                if ((size_t)cov >= max_cover) { // resize buffer
                    size_t new_cap = max_cover;
                    while ((size_t)cov >= new_cap) {
                        new_cap *= 2;
                    }
                    dist = fc_recalloc(dist, max_cover, new_cap, sizeof(int64_t), "coverage distribution");
                    max_cover = new_cap;
                }
                dist[(size_t)cov] += 1;
            }

            // summary stats including sparse (user-defined) coverage distribution
            _write_summary(NULL, wr->fh_summary, reg, stats,
                dist, max_cover, wr->thresholds, wr->n_thresholds);

            // update total stats - these get written on close
            wr->stats[0] = min(stats[0], wr->stats[0]);
            wr->stats[1] = max(stats[1], wr->stats[1]);
            wr->stats[2] += stats[2];
            wr->stats[3] += stats[3];
            // ..and total dist
            if (max_cover > wr->max_cover) {
                wr->dist = fc_recalloc(wr->dist, wr->max_cover, max_cover, sizeof(int64_t), "total coverage distribution");
                wr->max_cover = max_cover;
            }
            for (size_t i = 0; i < max_cover; ++i) {
                wr->dist[i] += dist[i];
            }
            free(dist);

            // process coverage into piecewise constant segments for output
            if (wr->per_base) {
                int64_t pos = reg->start;
                // cumulative coverage
                uint32_t cov_fwd = (uint32_t)w->diff_fwd[pos];
                uint32_t cov_rev = (uint32_t)w->diff_rev[pos];
                uint32_t cov_tot = cov_fwd + cov_rev;
                // current segment starts
                int64_t seg_start_fwd = reg->start;
                int64_t seg_start_rev = reg->start;
                int64_t seg_start_tot = reg->start;

                for (pos = reg->start + 1; pos < reg->end; ++pos) {
                    const uint32_t prev_f = cov_fwd;
                    const uint32_t prev_r = cov_rev;
                    const uint32_t prev_tot = prev_f + prev_r;

                    cov_fwd = (uint32_t)w->diff_fwd[pos];
                    cov_rev = (uint32_t)w->diff_rev[pos];
                    cov_tot = cov_fwd + cov_rev;

                    if ((cov_fwd != prev_f)) { // forward strand
                        _write_bedgraph_line(wr->fh_fwd, w->chrom, seg_start_fwd, pos, prev_f);
                        seg_start_fwd = pos;
                    }
                    if (cov_rev != prev_r) { // reverse strand
                        _write_bedgraph_line(wr->fh_rev, w->chrom, seg_start_rev, pos, prev_r);
                        seg_start_rev = pos;
                    }
                    if (cov_tot != prev_tot) { // total coverage
                        _write_bedgraph_line(wr->fh_tot, w->chrom, seg_start_tot, pos, prev_tot);
                        seg_start_tot = pos;
                    }
                }

                // close trailing segments
                if (seg_start_fwd < reg->end) {
                    _write_bedgraph_line(wr->fh_fwd, w->chrom, seg_start_fwd, reg->end, cov_fwd);
                }
                if (seg_start_rev < reg->end) {
                    _write_bedgraph_line(wr->fh_rev, w->chrom, seg_start_rev, reg->end, cov_rev);
                }
                if (seg_start_tot < reg->end) {
                    _write_bedgraph_line(wr->fh_tot, w->chrom, seg_start_tot, reg->end, cov_tot);
                }
            }

            if (!wr->whole_chrom) {
                free(region_name);
            }
            wr->cur_region = j + 1;

        } // for each region in BED file
    } // for each BED file
}

static void _reset_contig(cov_writer w, const int32_t tid)
{
    if (tid >= 0) {
        if (w->hdr == NULL || tid >= w->hdr->n_targets) {
            FC_LOG_FATAL("Invalid contig id %d (header has %d targets).", tid, w->hdr ? w->hdr->n_targets : 0);
        }
        w->tid = tid;
        w->contig_len = w->hdr->target_len[tid];
        if (w->contig_len < 0) {
            FC_LOG_FATAL("Invalid contig length %" PRId64 " for %s.", w->contig_len, w->hdr->target_name[tid]);
        }
        if ((uint64_t)w->contig_len > ((uint64_t)SIZE_MAX - 1u)) {
            FC_LOG_FATAL("Contig %s length %" PRId64 " cannot be represented on this platform.",
                w->hdr->target_name[tid], w->contig_len);
        }
        w->chrom = w->hdr->target_name[tid];

        size_t required_buf = (size_t)w->contig_len + 1u;
        if (required_buf > (SIZE_MAX / sizeof(int32_t))) {
            FC_LOG_FATAL("Coverage buffer for contig %s would overflow addressable memory.",
                w->hdr->target_name[tid]);
        }

        if (required_buf > w->buf_size) {
            // grow buffer and clear everything
            w->buf_size = required_buf;
            w->diff_fwd = fc_realloc(w->diff_fwd, w->buf_size * sizeof(int32_t), "coverage buffer");
            w->diff_rev = fc_realloc(w->diff_rev, w->buf_size * sizeof(int32_t), "coverage buffer");
            memset(w->diff_fwd, 0, w->buf_size * sizeof(int32_t));
            memset(w->diff_rev, 0, w->buf_size * sizeof(int32_t));
        } else {
            // buffer is good size, zero what we need
            memset(w->diff_fwd, 0, required_buf * sizeof(int32_t));
            memset(w->diff_rev, 0, required_buf * sizeof(int32_t));
        }
    } else {
        w->tid = -1;
        w->contig_len = 0;
        w->chrom = NULL;
    }
}

cov_writer_region init_coverage_writer_region(
    const char *out_dir, const char *name, bool per_base, bool by_strand, const char *bed_fname,
    uint32_t *thresholds, size_t n_thresholds, uint32_t segment_length,
    const bam_hdr_t *hdr, const htsThreadPool *pool)
{

    if (mkdir_hier((char *)out_dir, true) == -1) {
        FC_LOG_FATAL("Failed to create coverage output directory '%s'.", out_dir);
    }

    _cov_writer_region *w = fc_alloc(1, sizeof(_cov_writer_region), "coverage writer region");
    w->name = fc_strdup(name, "coverage region name");
    w->per_base = per_base;
    if (bed_fname != NULL) {
        w->whole_chrom = false;
        w->bed = init_bed(bed_fname, hdr);
    } else {
        // no BED file, use full genome
        if (segment_length != 0) {
            w->whole_chrom = false;
            w->bed = init_bed_from_sam(hdr, segment_length);
        } else {
            w->whole_chrom = true;
            w->bed = init_bed_from_sam(hdr, 0);
        }
    }
    w->cur_region = 0;

    // sparse CDF thresholds
    if (n_thresholds > 0 && thresholds) {
        w->thresholds = fc_alloc(n_thresholds, sizeof *w->thresholds, "thresholds");
        memcpy(w->thresholds, thresholds, n_thresholds * sizeof *w->thresholds);
        w->n_thresholds = n_thresholds;
        qsort(w->thresholds, w->n_thresholds, sizeof *w->thresholds, cmp_u32);
    } else {
        // use default
        w->n_thresholds = 6;
        w->thresholds = fc_alloc(w->n_thresholds, sizeof *w->thresholds, "default thresholds");
        uint32_t tmp[] = { 1, 5, 10, 20, 30, 40 };
        for (size_t i = 0; i < w->n_thresholds; ++i) {
            w->thresholds[i] = tmp[i];
        }
    }
    // min, max, total bases, num. positions
    w->stats[0] = INT64_MAX;
    w->stats[1] = 0;
    w->stats[2] = 0;
    w->stats[3] = 0;
    // dist needs to be large enough to hold all thresholds, start with at least 256 to avoid reallocating
    w->max_cover = max(1u << 8, w->thresholds ? 1 + w->thresholds[w->n_thresholds - 1] : 0);
    w->dist = fc_alloc(w->max_cover, sizeof(int64_t), "coverage distribution");

    // per_base outputs are the main event, the full coverage traces.
    if (w->per_base) {
        const char *suffixes[3] = { ".fwd.bed.gz", ".rev.bed.gz", ".bed.gz" };
        char **fn[3] = { &w->fn_fwd, &w->fn_rev, &w->fn_tot };
        BGZF **fh[3] = { &w->fh_fwd, &w->fh_rev, &w->fh_tot };

        int k = by_strand ? 0 : 2; // only need last is not by_strand
        for (int i = k; i < 3; ++i) {
            *fn[i] = fc_asprintf("bedgraph output path", "%s/%s%s", out_dir, name, suffixes[i]);
            *fh[i] = bgzf_open(*fn[i], "w1");
            if (NULL == *fh[i]) {
                FC_LOG_FATAL("Cannot open bedgraph output '%s'", *fn[i]);
            }
            bgzf_index_build_init(*fh[i]);
            if (NULL != pool && NULL != pool->pool) {
                bgzf_thread_pool(*fh[i], pool->pool, 0);
            }
        }
    }

    // summary file {name}.mosdepth.summary.txt
    {
        char *suffix = ".summary.txt";
        char *fname = fc_asprintf("summary output path", "%s/%s%s", out_dir, name, suffix);
        w->fh_summary = fopen(fname, "w");
        if (NULL == w->fh_summary) {
            FC_LOG_FATAL("Cannot open summary output '%s'", fname);
        }
        // TODO: mosdepth allows mean to be changed to median
        fprintf(w->fh_summary, "chrom\tstart\tend\tlength\tbases\tmean\tmin\tmax");
        for (size_t i = 0; i < w->n_thresholds; ++i) {
            fprintf(w->fh_summary, "\t%dx", w->thresholds[i]);
        }
        fprintf(w->fh_summary, "\n");
        free(fname);
    }

    // distribution of coverage {name}.mosdepth.dist.txt
    {
        char *suffix = ".dist.txt";
        char *fname = fc_asprintf("distribution output path", "%s/%s%s", out_dir, name, suffix);
        w->fh_dist = fopen(fname, "w");
        if (NULL == w->fh_dist) {
            FC_LOG_FATAL("Cannot open summary output '%s'", fname);
        }
        free(fname);
    }

    // NOTE: mosdepth can emit a separate ".thresholds.txt" file, but fastcat does not.
    // Threshold CDF values are already written into ".summary.txt" as wide "Nx" columns,
    // and per-threshold counts are derived from the same cumulative distribution.

    return w;
}

void destroy_coverage_writer_region(cov_writer_region w)
{
    if (NULL == w) {
        return;
    }

    // data structures associated with BED regions
    if (NULL != w->bed) {
        destroy_bed(w->bed);
    }

    // distribution and thresholds files
    // make a region corresponding to the whole BED for writing distribution
    _bed_region reg = { w->name, 0, 0, w->stats[3] };
    _write_summary(w->fh_dist, w->fh_summary, &reg, w->stats,
        w->dist, w->max_cover, w->thresholds, w->n_thresholds);
    fclose(w->fh_dist);
    // fclose(w->fh_thresh);
    fclose(w->fh_summary);

    // create BED index files if we wrote per-base coverage
    if (w->per_base) {
        BGZF *fhs[3] = { w->fh_fwd, w->fh_rev, w->fh_tot };
        char **fns[3] = { &w->fn_fwd, &w->fn_rev, &w->fn_tot };
        for (int i = 0; i < 3; ++i) {
            if (fhs[i] == NULL) {
                continue; // from by_strand == false
            }
            if (bgzf_index_dump(fhs[i], *fns[i], ".csi") < 0) {
                FC_LOG_FATAL("Cannot write index for '%s'", *fns[i]);
            }
            if (bgzf_close(fhs[i]) < 0) {
                FC_LOG_FATAL("Cannot close '%s'", *fns[i]);
            }
            free(*fns[i]);
        }
    }

    free(w->name);
    free(w->thresholds);
    free(w->dist);

    free(w);
}

cov_writer init_coverage_writer(
    const char *out_dir, bool per_base, bool use_cigar,
    int exclude_flags, int include_flags, bool by_strand,
    const bam_hdr_t *hdr, const htsThreadPool *pool,
    const char **bed_files, const char **bed_names, size_t n_beds,
    uint32_t *thresholds, size_t n_thresholds,
    uint32_t *segments, size_t n_segments)
{

    // check outputs are writeable, we don't use mkdir_hier here because otherwise
    // the global coverage writer would fail to make the TLD.
    if (can_make_dir(out_dir) < 1) {
        FC_LOG_FATAL(
            "Cannot create output directory '%s'. Check location is writeable and directory does not exist.",
            out_dir);
    }

    kstring_t ks = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd((bam_hdr_t *)hdr, "SO", &ks) >= 0) {
        if (strcmp(ks.s, "coordinate") != 0) {
            FC_LOG_FATAL("BAM sort order is '%s', must be 'coordinate'. Cannot calculate coverage data.", ks.s);
        }
    } else {
        FC_LOG_FATAL("BAM header missing 'SO' (sort order) tag. Cannot calculate coverage data.");
    }
    free(ks.s);

    cov_writer w = (cov_writer)fc_alloc(1, sizeof(_cov_writer), "coverage writer");

    w->use_cigar = use_cigar;
    w->exclude_flags = exclude_flags == -1 ? 1796 : exclude_flags; // default from mosdepth: unmapped, not primary, failed QC, duplicate
    w->include_flags = include_flags == -1 ? 0 : include_flags; // default from mosdepth: no flags
    w->buf_size = 0;
    w->diff_fwd = NULL;
    w->diff_rev = NULL;
    w->hdr = hdr;
    _reset_contig(w, -1);

    w->n_beds = 1 + n_segments + n_beds;
    w->writers = fc_alloc(w->n_beds, sizeof(cov_writer_region), "coverage writers");
    // first entry is complete genome
    w->writers[0] = init_coverage_writer_region(
        out_dir, "global", per_base, by_strand, NULL,
        thresholds, n_thresholds, 0,
        hdr, pool);
    // fixed length segments
    for (size_t i = 0; i < n_segments; ++i) {
        char seg_name[32]; // largest uint32 is 10 digits
        FC_SNPRINTF_CHECK(seg_name, sizeof(seg_name), "segments_%u", segments[i]);
        char *seg_out_dir = fc_asprintf("segment output directory path", "%s/segments_%u", out_dir, segments[i]);
        w->writers[1 + i] = init_coverage_writer_region(
            seg_out_dir, seg_name, per_base, by_strand, NULL,
            thresholds, n_thresholds, segments[i],
            hdr, pool);
        free(seg_out_dir);
    }
    // followed by the BED files
    for (size_t i = 0; i < n_beds; ++i) {
        if (NULL != bed_files[i]) {
            const char *bed_name = bed_names[i];
            char *bed_out_dir = fc_asprintf("bed output directory path", "%s/%s", out_dir, bed_name);
            w->writers[1 + n_segments + i] = init_coverage_writer_region(
                bed_out_dir, bed_name, per_base, by_strand, bed_files[i],
                thresholds, n_thresholds, 0,
                hdr, pool);
            free(bed_out_dir);
        }
    }
    return w;
}

void destroy_coverage_writer(cov_writer w)
{
    if (NULL == w) {
        return;
    }

    // emit remaining coverage
    if (w->tid >= 0) {
        _flush_contig(w);
    }

    // flush any tail-end regions in BED files
    w->tid = w->hdr->n_targets; // all BED regions will have a tid < n_targets
    _fill_skipped_regions(w);

    for (size_t i = 0; i < w->n_beds; ++i) {
        if (w->writers[i] == NULL) {
            continue;
        }
        destroy_coverage_writer_region(w->writers[i]);
    }
    free(w->writers);

    free(w->diff_fwd);
    free(w->diff_rev);
    free(w);
}

void coverage_process(cov_writer w, const bam1_t *b)
{
    if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) {
        return;
    }
    // apply exclude flags THEN include flags (matching mosdepth)
    if (b->core.flag & w->exclude_flags) {
        return;
    }
    if (w->include_flags && !(b->core.flag & w->include_flags)) {
        return;
    }

    // flush previous contig
    if (w->tid != b->core.tid) {
        _flush_contig(w);
        _reset_contig(w, b->core.tid);
    }

    // select array to update
    const int is_rev = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
    int32_t *diff = is_rev ? w->diff_rev : w->diff_fwd;

    // either use rstart/rend or studiously examine cigar
    if (!w->use_cigar) {
        int64_t rstart = (int64_t)b->core.pos;
        int64_t rend = (int64_t)bam_endpos(b);
        if (rstart < 0 || rend < rstart || rstart >= w->contig_len || rend > w->contig_len) {
            coverage_fail_invalid_alignment(w, b, "reference span falls outside contig bounds");
        }
        diff[rstart] += 1;
        diff[rend] -= 1;
    } else {
        const uint32_t *cigar = bam_get_cigar(b);
        int64_t rstart = (int64_t)b->core.pos;
        if (rstart < 0 || rstart > w->contig_len) {
            coverage_fail_invalid_alignment(w, b, "reference start falls outside contig bounds");
        }
        for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
            const uint32_t op = cigar[i];
            int64_t len = (int64_t)bam_cigar_oplen(op);
            if (is_ref_and_query(op)) {
                if (len < 0 || rstart + len > w->contig_len) {
                    coverage_fail_invalid_alignment(w, b, "CIGAR consumes reference past contig end");
                }
                diff[rstart] += 1;
                diff[rstart + len] -= 1;
                rstart += len;
            } else if (consumes_ref_only(op)) {
                if (len < 0 || rstart + len > w->contig_len) {
                    coverage_fail_invalid_alignment(w, b, "CIGAR skips past contig end");
                }
                rstart += len; // advance reference, no coverage
            } else {
                // do not advance reference
            }
        }
    }
}
