#include "core/fc_report.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "htslib/kstring.h"

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "core/fc_report_template_data.h"
#include "sink/common.h"
#include "version.h"

#define FC_REPORT_LENGTH_BINS 100
#define FC_REPORT_QUALITY_MIN 4.0
#define FC_REPORT_QUALITY_MAX 30.0
#define FC_REPORT_QUALITY_WIDTH 0.2
#define FC_REPORT_QUALITY_BINS 130
#define FC_REPORT_YIELD_THINNING 100
#define FC_REPORT_LQ2D_COLOR_QUANTILE 0.995
#define FC_REPORT_NAME "fastcat_report.html"

typedef struct {
    double start;
    double end;
    double count;
} report_hist_bin;

typedef struct {
    report_hist_bin *bins;
    size_t n;
    size_t cap;
    double total_count;
} report_hist;

static void report_kputs(kstring_t *ks, const char *text)
{
    if (text && kputs(text, ks) < 0) {
        FC_LOG_FATAL("Failed to append report string.");
    }
}

static void report_kputsn(kstring_t *ks, const char *text, size_t len)
{
    if (len > 0 && kputsn(text, len, ks) < 0) {
        FC_LOG_FATAL("Failed to append report string.");
    }
}

static void report_kputc(kstring_t *ks, int c)
{
    if (kputc(c, ks) < 0) {
        FC_LOG_FATAL("Failed to append report string.");
    }
}

static void report_ksprintf(kstring_t *ks, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int rc = kvsprintf(ks, fmt, ap);
    va_end(ap);
    if (rc < 0) {
        FC_LOG_FATAL("Failed to format report string.");
    }
}

static char *report_steal(kstring_t *ks)
{
    if (!ks->s) {
        report_kputs(ks, "");
    }
    char *data = ks->s;
    ks->s = NULL;
    ks->l = 0;
    ks->m = 0;
    return data;
}

static void report_hist_destroy(report_hist *hist)
{
    free(hist->bins);
    hist->bins = NULL;
    hist->n = 0;
    hist->cap = 0;
    hist->total_count = 0.0;
}

static bool report_hist_push(report_hist *hist, double start, double end, double count)
{
    if (!isfinite(start) || !isfinite(end) || !isfinite(count) || count < 0.0) {
        return false;
    }
    if (count == 0.0) {
        return true;
    }
    if (hist->n == hist->cap) {
        size_t new_cap = hist->cap ? hist->cap * 2 : 256;
        hist->bins = fc_realloc(hist->bins, new_cap * sizeof(*hist->bins), "report histogram bins");
        hist->cap = new_cap;
    }
    hist->bins[hist->n++] = (report_hist_bin) {
        .start = start,
        .end = end,
        .count = count
    };
    hist->total_count += count;
    return true;
}

static bool report_hist_from_read_stats(const read_stats *stats, report_hist *hist)
{
    if (!stats || !hist) {
        return false;
    }

    for (size_t i = 0; i < stats->n; ++i) {
        size_t count = stats->counts[i];
        if (count == 0) {
            continue;
        }

        double start = 0.0;
        double end = 0.0;
        if (stats->width == 0.0f) {
            start = (double)stats->edges[i];
            end = (double)stats->edges[i + 1];
        } else {
            start = (double)i * (double)stats->width;
            end = (double)(i + 1) * (double)stats->width;
        }
        if (!report_hist_push(hist, start, end, (double)count)) {
            return false;
        }
    }
    return true;
}

static double read_stats_total_count(const read_stats *stats)
{
    if (!stats) {
        return 0.0;
    }
    double total = 0.0;
    for (size_t i = 0; i < stats->n; ++i) {
        total += (double)stats->counts[i];
    }
    return total;
}

static void select_report_stats(
    const fc_stats_engine *engine,
    bool unmapped_override,
    const read_stats **length,
    const read_stats **quality)
{
    *length = NULL;
    *quality = NULL;
    if (!engine) {
        return;
    }

    if (unmapped_override) {
        *length = engine->length_stats_unmapped;
        *quality = engine->qual_stats_unmapped;
    } else {
        *length = engine->length_stats;
        *quality = engine->qual_stats;
    }
}

static void append_json_string(kstring_t *ks, const char *text)
{
    report_kputs(ks, "\"");
    if (text) {
        for (const unsigned char *p = (const unsigned char *)text; *p; ++p) {
            switch (*p) {
            case '"':
                report_kputs(ks, "\\\"");
                break;
            case '\\':
                report_kputs(ks, "\\\\");
                break;
            case '\b':
                report_kputs(ks, "\\b");
                break;
            case '\f':
                report_kputs(ks, "\\f");
                break;
            case '\n':
                report_kputs(ks, "\\n");
                break;
            case '\r':
                report_kputs(ks, "\\r");
                break;
            case '\t':
                report_kputs(ks, "\\t");
                break;
            default:
                if (*p < 0x20) {
                    report_ksprintf(ks, "\\u%04x", *p);
                } else {
                    report_kputc(ks, *p);
                }
                break;
            }
        }
    }
    report_kputs(ks, "\"");
}

static void append_html_escaped(kstring_t *ks, const char *text)
{
    if (!text) {
        return;
    }
    for (const char *p = text; *p; ++p) {
        switch (*p) {
        case '&':
            report_kputs(ks, "&amp;");
            break;
        case '<':
            report_kputs(ks, "&lt;");
            break;
        case '>':
            report_kputs(ks, "&gt;");
            break;
        case '"':
            report_kputs(ks, "&quot;");
            break;
        case '\'':
            report_kputs(ks, "&#39;");
            break;
        default:
            report_kputc(ks, *p);
            break;
        }
    }
}

static char *html_escape(const char *text)
{
    kstring_t ks = { 0, 0, NULL };
    append_html_escaped(&ks, text);
    return report_steal(&ks);
}

static char *build_run_config_rows_html(const fc_report_config_entry *entries, size_t count)
{
    kstring_t ks = { 0, 0, NULL };
    if (!entries || count == 0) {
        report_kputs(&ks, "<tr><td>run configuration</td><td>not available</td></tr>\n");
        return report_steal(&ks);
    }

    for (size_t i = 0; i < count; ++i) {
        report_kputs(&ks, "<tr><td>");
        append_html_escaped(&ks, entries[i].key ? entries[i].key : "(not set)");
        report_kputs(&ks, "</td><td>");
        append_html_escaped(&ks, entries[i].value ? entries[i].value : "(not set)");
        report_kputs(&ks, "</td></tr>\n");
    }
    return report_steal(&ks);
}

static void replace_all(char **html, const char *needle, const char *replacement)
{
    kstring_t ks = { 0, 0, NULL };
    const char *cursor = *html;
    size_t needle_len = strlen(needle);
    while (true) {
        const char *hit = strstr(cursor, needle);
        if (!hit) {
            report_kputs(&ks, cursor);
            break;
        }
        report_kputsn(&ks, cursor, (size_t)(hit - cursor));
        report_kputs(&ks, replacement ? replacement : "");
        cursor = hit + needle_len;
    }
    free(*html);
    *html = report_steal(&ks);
}

static char *inflate_report_template(void)
{
    if (fc_report_template_html_len > SIZE_MAX - 1 || fc_report_template_html_len > UINT_MAX
        || fc_report_template_html_gz_len > UINT_MAX) {
        FC_LOG_ERROR("Embedded report template is too large to inflate.");
        return NULL;
    }

    char *html = fc_alloc(fc_report_template_html_len + 1, sizeof(*html), "report template");
    z_stream stream = { 0 };
    stream.next_in = (Bytef *)fc_report_template_html_gz;
    stream.avail_in = (uInt)fc_report_template_html_gz_len;
    stream.next_out = (Bytef *)html;
    stream.avail_out = (uInt)fc_report_template_html_len;

    int zret = inflateInit2(&stream, 16 + MAX_WBITS);
    if (zret != Z_OK) {
        FC_LOG_ERROR("Failed to initialise report template inflater.");
        free(html);
        return NULL;
    }
    zret = inflate(&stream, Z_FINISH);
    int zend = inflateEnd(&stream);

    if (zret != Z_STREAM_END || zend != Z_OK || stream.total_out != fc_report_template_html_len) {
        FC_LOG_ERROR("Failed to inflate embedded report template.");
        free(html);
        return NULL;
    }

    html[fc_report_template_html_len] = '\0';
    return html;
}

static double hist_mean_midpoint(const report_hist *hist)
{
    if (!hist || hist->total_count <= 0.0) {
        return 0.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < hist->n; ++i) {
        double mid = 0.5 * (hist->bins[i].start + hist->bins[i].end);
        sum += mid * hist->bins[i].count;
    }
    return sum / hist->total_count;
}

static double hist_mean_start(const report_hist *hist)
{
    if (!hist || hist->total_count <= 0.0) {
        return 0.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < hist->n; ++i) {
        sum += hist->bins[i].start * hist->bins[i].count;
    }
    return sum / hist->total_count;
}

static double hist_median_start(const report_hist *hist)
{
    if (!hist || hist->total_count <= 0.0 || hist->n == 0) {
        return 0.0;
    }
    double target = 0.5 * hist->total_count;
    double cumsum = 0.0;
    for (size_t i = 0; i < hist->n; ++i) {
        cumsum += hist->bins[i].count;
        if (cumsum >= target) {
            return hist->bins[i].start;
        }
    }
    return hist->bins[hist->n - 1].start;
}

static void append_number_array(kstring_t *ks, const double *values, size_t n)
{
    report_kputs(ks, "[");
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) {
            report_kputs(ks, ",");
        }
        double value = isfinite(values[i]) ? values[i] : 0.0;
        report_ksprintf(ks, "%.17g", value);
    }
    report_kputs(ks, "]");
}

static void append_base64(kstring_t *ks, const uint8_t *data, size_t n)
{
    static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    for (size_t i = 0; i < n; i += 3) {
        uint32_t a = data[i];
        uint32_t b = (i + 1 < n) ? data[i + 1] : 0;
        uint32_t c = (i + 2 < n) ? data[i + 2] : 0;
        uint32_t triple = (a << 16) | (b << 8) | c;

        report_kputc(ks, table[(triple >> 18) & 0x3f]);
        report_kputc(ks, table[(triple >> 12) & 0x3f]);
        report_kputc(ks, (i + 1 < n) ? table[(triple >> 6) & 0x3f] : '=');
        report_kputc(ks, (i + 2 < n) ? table[triple & 0x3f] : '=');
    }
}

static void append_lq2d_empty_plot_json(kstring_t *json)
{
    report_kputs(json, "\"length_qscore\":{");
    report_kputs(json, "\"image\":null,\"x\":0,\"y\":0,\"dw\":1,\"dh\":1,");
    report_kputs(json, "\"color_high\":1,\"subtitle\":");
    append_json_string(json, "No length/qscore data available");
    report_kputs(json, "}");
}

static void append_lq2d_image_ndarray_json(
    kstring_t *json,
    const fc_length_qscore_2d_matrix *matrix)
{
    if (matrix->length_bins != 0 && matrix->qscore_bins > SIZE_MAX / matrix->length_bins) {
        FC_LOG_FATAL("Length/qscore report image is too large.");
    }
    size_t n = matrix->length_bins * matrix->qscore_bins;
    if (n > SIZE_MAX / sizeof(double)) {
        FC_LOG_FATAL("Length/qscore report image is too large.");
    }

    size_t byte_count = n * sizeof(double);
    uint8_t *bytes = fc_alloc(byte_count, sizeof(*bytes), "report length/qscore image bytes");
    size_t out = 0;
    for (size_t len_bin = 0; len_bin < matrix->length_bins; len_bin++) {
        for (size_t q = 0; q < matrix->qscore_bins; q++) {
            double count = matrix->counts[q * matrix->length_bins + len_bin];
            if (!isfinite(count) || count < 0.0) {
                count = 0.0;
            }

            uint64_t bits = 0;
            memcpy(&bits, &count, sizeof(bits));
            for (size_t byte = 0; byte < sizeof(bits); byte++) {
                bytes[out++] = (uint8_t)((bits >> (8 * byte)) & 0xff);
            }
        }
    }

    report_kputs(json, "{\"type\":\"ndarray\",\"array\":{\"type\":\"bytes\",\"data\":\"");
    append_base64(json, bytes, byte_count);
    report_kputs(json, "\"},\"shape\":[");
    report_ksprintf(json, "%zu,%zu", matrix->length_bins, matrix->qscore_bins);
    report_kputs(json, "],\"dtype\":\"float64\",\"order\":\"little\"}");
    free(bytes);
}

static double lq2d_color_high_quantile(const fc_length_qscore_2d_matrix *matrix)
{
    if (matrix->length_bins != 0 && matrix->qscore_bins > SIZE_MAX / matrix->length_bins) {
        FC_LOG_FATAL("Length/qscore report image is too large.");
    }
    size_t n = matrix->length_bins * matrix->qscore_bins;
    if (n == 0) {
        return 1.0;
    }

    double *values = fc_alloc(n, sizeof(*values), "report length/qscore color range");
    for (size_t i = 0; i < n; i++) {
        double value = matrix->counts[i];
        values[i] = (isfinite(value) && value > 0.0) ? value : 0.0;
    }

    qsort(values, n, sizeof(*values), cmp_double);
    size_t idx = (size_t)llround((double)(n - 1) * FC_REPORT_LQ2D_COLOR_QUANTILE);
    if (idx >= n) {
        idx = n - 1;
    }

    double high = values[idx];
    if (!isfinite(high) || high <= 0.0) {
        high = values[n - 1];
    }
    free(values);

    return (isfinite(high) && high > 0.0) ? high : 1.0;
}

static void build_length_qscore_plot_json(
    kstring_t *json,
    const fc_length_qscore_2d *model,
    size_t barcode)
{
    if (!model || !fc_length_qscore_2d_active(model, barcode)) {
        append_lq2d_empty_plot_json(json);
        return;
    }

    fc_length_qscore_2d_matrix matrix = { 0 };
    if (fc_length_qscore_2d_materialize(model, barcode, &matrix) != 0 || !matrix.counts) {
        fc_length_qscore_2d_matrix_destroy(&matrix);
        append_lq2d_empty_plot_json(json);
        return;
    }

    double total = 0.0;
    size_t n = matrix.length_bins * matrix.qscore_bins;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(matrix.counts[i]) && matrix.counts[i] > 0.0) {
            total += matrix.counts[i];
        }
    }

    char subtitle[160];
    FC_SNPRINTF_CHECK(subtitle, sizeof(subtitle), "Estimated counts, n=%.0f, length cap %.3gkb",
        total, matrix.length_upper / 1000.0);

    double color_high = lq2d_color_high_quantile(&matrix);
    report_kputs(json, "\"length_qscore\":{");
    report_kputs(json, "\"image\":");
    append_lq2d_image_ndarray_json(json, &matrix);
    report_kputs(json, ",\"x\":0,\"y\":0,\"dw\":");
    report_ksprintf(json, "%.17g", matrix.qscore_width * (double)matrix.qscore_bins);
    report_kputs(json, ",\"dh\":");
    report_ksprintf(json, "%.17g", matrix.length_upper / 1000.0);
    report_kputs(json, ",\"color_high\":");
    report_ksprintf(json, "%.17g", color_high > 0.0 ? color_high : 1.0);
    report_kputs(json, ",\"subtitle\":");
    append_json_string(json, subtitle);
    report_kputs(json, "}");

    fc_length_qscore_2d_matrix_destroy(&matrix);
}

static void build_quality_plot_json(kstring_t *json, const report_hist *quality)
{
    double top[FC_REPORT_QUALITY_BINS] = { 0 };
    double left[FC_REPORT_QUALITY_BINS] = { 0 };
    double right[FC_REPORT_QUALITY_BINS] = { 0 };

    for (size_t i = 0; i < FC_REPORT_QUALITY_BINS; ++i) {
        left[i] = FC_REPORT_QUALITY_MIN + (double)i * FC_REPORT_QUALITY_WIDTH;
        right[i] = left[i] + FC_REPORT_QUALITY_WIDTH;
    }
    for (size_t i = 0; i < quality->n; ++i) {
        double value = quality->bins[i].start;
        if (value < FC_REPORT_QUALITY_MIN || value > FC_REPORT_QUALITY_MAX) {
            continue;
        }
        size_t bin = (value == FC_REPORT_QUALITY_MAX)
            ? FC_REPORT_QUALITY_BINS - 1
            : (size_t)floor((value - FC_REPORT_QUALITY_MIN) / FC_REPORT_QUALITY_WIDTH);
        if (bin < FC_REPORT_QUALITY_BINS) {
            top[bin] += quality->bins[i].count;
        }
    }

    char subtitle[128];
    FC_SNPRINTF_CHECK(subtitle, sizeof(subtitle), "Mean: %.1f. Median: %.1f",
        hist_mean_midpoint(quality), hist_median_start(quality));

    report_kputs(json, "\"quality\":{");
    report_kputs(json, "\"top\":");
    append_number_array(json, top, FC_REPORT_QUALITY_BINS);
    report_kputs(json, ",\"left\":");
    append_number_array(json, left, FC_REPORT_QUALITY_BINS);
    report_kputs(json, ",\"right\":");
    append_number_array(json, right, FC_REPORT_QUALITY_BINS);
    report_kputs(json, ",\"subtitle\":");
    append_json_string(json, subtitle);
    report_kputs(json, "}");
}

static void build_length_plot_json(kstring_t *json, const report_hist *length)
{
    double top[FC_REPORT_LENGTH_BINS] = { 0 };
    double left[FC_REPORT_LENGTH_BINS] = { 0 };
    double right[FC_REPORT_LENGTH_BINS] = { 0 };

    double min_kb = (length->n > 0) ? length->bins[0].start / 1000.0 : 0.0;
    double max_kb = (length->n > 0) ? length->bins[length->n - 1].start / 1000.0 : 1.0;
    if (max_kb <= min_kb) {
        max_kb = min_kb + 1.0;
    }
    double width = (max_kb - min_kb) / (double)FC_REPORT_LENGTH_BINS;
    for (size_t i = 0; i < FC_REPORT_LENGTH_BINS; ++i) {
        left[i] = min_kb + (double)i * width;
        right[i] = left[i] + width;
    }
    for (size_t i = 0; i < length->n; ++i) {
        double value = length->bins[i].start / 1000.0;
        size_t bin = (value >= max_kb)
            ? FC_REPORT_LENGTH_BINS - 1
            : (size_t)floor((value - min_kb) / width);
        if (bin < FC_REPORT_LENGTH_BINS) {
            top[bin] += length->bins[i].count;
        }
    }

    char subtitle[192];
    uint64_t min_len = (length->n > 0) ? (uint64_t)llround(length->bins[0].start) : 0;
    uint64_t max_len = (length->n > 0) ? (uint64_t)llround(length->bins[length->n - 1].start) : 0;
    FC_SNPRINTF_CHECK(subtitle, sizeof(subtitle),
        "Mean: %.1f. Median: %.1f. Min: %" PRIu64 ". Max: %" PRIu64,
        hist_mean_start(length), hist_median_start(length), min_len, max_len);

    report_kputs(json, "\"length\":{");
    report_kputs(json, "\"top\":");
    append_number_array(json, top, FC_REPORT_LENGTH_BINS);
    report_kputs(json, ",\"left\":");
    append_number_array(json, left, FC_REPORT_LENGTH_BINS);
    report_kputs(json, ",\"right\":");
    append_number_array(json, right, FC_REPORT_LENGTH_BINS);
    report_kputs(json, ",\"subtitle\":");
    append_json_string(json, subtitle);
    report_kputs(json, "}");
}

static double compute_n50_kb(const report_hist *length, double total_yield)
{
    if (!length || length->n == 0 || total_yield <= 0.0) {
        return 0.0;
    }
    double target = total_yield * 0.5;
    double acc = 0.0;
    for (size_t i = length->n; i > 0; --i) {
        const report_hist_bin *bin = &length->bins[i - 1];
        acc += bin->start * bin->count;
        if (acc >= target) {
            return bin->start / 1000.0;
        }
    }
    return length->bins[0].start / 1000.0;
}

static void build_yield_plot_json(kstring_t *json, const report_hist *length)
{
    double *x = NULL;
    double *y = NULL;
    size_t out_n = 0;
    size_t out_cap = 0;
    double *cumsum = NULL;
    double total_yield = 0.0;

    if (length->n > 0) {
        cumsum = fc_alloc(length->n, sizeof(*cumsum), "report yield cumsum");
        for (size_t i = length->n; i > 0; --i) {
            const report_hist_bin *bin = &length->bins[i - 1];
            total_yield += bin->start * bin->count;
            cumsum[i - 1] = total_yield;
        }

        size_t step = (length->n > FC_REPORT_YIELD_THINNING)
            ? length->n / FC_REPORT_YIELD_THINNING
            : 1;
        for (size_t i = 0; i < length->n; i += step) {
            if (out_n == out_cap) {
                out_cap = out_cap ? out_cap * 2 : 128;
                x = fc_realloc(x, out_cap * sizeof(*x), "report yield x");
                y = fc_realloc(y, out_cap * sizeof(*y), "report yield y");
            }
            x[out_n] = length->bins[i].start / 1000.0;
            y[out_n] = cumsum[i] / 1e9;
            out_n++;
        }
        if (out_n == 0 || x[out_n - 1] != length->bins[length->n - 1].start / 1000.0) {
            if (out_n == out_cap) {
                out_cap = out_cap ? out_cap * 2 : 128;
                x = fc_realloc(x, out_cap * sizeof(*x), "report yield x");
                y = fc_realloc(y, out_cap * sizeof(*y), "report yield y");
            }
            x[out_n] = length->bins[length->n - 1].start / 1000.0;
            y[out_n] = cumsum[length->n - 1] / 1e9;
            out_n++;
        }
    }

    char subtitle[160];
    FC_SNPRINTF_CHECK(subtitle, sizeof(subtitle), "Total yield: %.3g Gb.  N50: %.3gkb",
        total_yield / 1e9, compute_n50_kb(length, total_yield));

    report_kputs(json, "\"yield\":{");
    report_kputs(json, "\"x\":");
    append_number_array(json, x, out_n);
    report_kputs(json, ",\"y\":");
    append_number_array(json, y, out_n);
    report_kputs(json, ",\"subtitle\":");
    append_json_string(json, subtitle);
    report_kputs(json, "}");

    free(cumsum);
    free(x);
    free(y);
}

static char *build_report_data_json(
    const report_hist *length,
    const report_hist *quality,
    const fc_length_qscore_2d *length_qscore,
    size_t barcode)
{
    kstring_t json = { 0, 0, NULL };
    report_kputs(&json, "{");
    build_quality_plot_json(&json, quality);
    report_kputs(&json, ",");
    build_length_plot_json(&json, length);
    report_kputs(&json, ",");
    build_yield_plot_json(&json, length);
    report_kputs(&json, ",");
    build_length_qscore_plot_json(&json, length_qscore, barcode);
    report_kputs(&json, "}");
    return report_steal(&json);
}

static char *render_report_html(
    const char *data_json,
    const char *sample,
    const char *sample_type,
    const char *barcode,
    const fc_report_config_entry *run_config,
    size_t run_config_count)
{
    char *html = inflate_report_template();
    if (!html) {
        return NULL;
    }

    char *sample_html = html_escape(sample ? sample : "fastcat");
    char *type_html = html_escape(sample_type ? sample_type : "fastcat");
    char *barcode_html = html_escape(barcode ? barcode : "unclassified");
    char *version_html = html_escape(argp_program_version ? argp_program_version : "unknown");
    char *run_config_rows_html = build_run_config_rows_html(run_config, run_config_count);

    replace_all(&html, "__FASTCAT_REPORT_DATA__", data_json);
    replace_all(&html, "__FASTCAT_SAMPLE_ALIAS_HTML__", sample_html);
    replace_all(&html, "__FASTCAT_SAMPLE_TYPE_HTML__", type_html);
    replace_all(&html, "__FASTCAT_SAMPLE_BARCODE_HTML__", barcode_html);
    replace_all(&html, "__FASTCAT_VERSION_HTML__", version_html);
    replace_all(&html, "__FASTCAT_RUN_CONFIG_ROWS__", run_config_rows_html);

    free(sample_html);
    free(type_html);
    free(barcode_html);
    free(version_html);
    free(run_config_rows_html);
    return html;
}

static int write_text_file(const char *path, const char *text)
{
    ensure_parent_dir_exists(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        FC_LOG_ERROR("Cannot open report '%s' for writing: %s", path, strerror(errno));
        return -1;
    }
    if (fputs(text, fp) == EOF) {
        FC_LOG_ERROR("Failed writing report '%s'.", path);
        fclose(fp);
        return -1;
    }
    if (fclose(fp) != 0) {
        FC_LOG_ERROR("Failed closing report '%s'.", path);
        return -1;
    }
    return 0;
}

static int write_report_for_output_dir(
    const char *output_dir,
    const char *sample,
    const char *sample_type,
    const char *barcode,
    const read_stats *length_stats,
    const read_stats *quality_stats,
    const fc_report_config_entry *run_config,
    size_t run_config_count,
    const fc_length_qscore_2d *length_qscore,
    size_t barcode_index)
{
    int rc = -1;
    report_hist length = { 0 };
    report_hist quality = { 0 };
    char *report_path = fc_join_path(output_dir, FC_REPORT_NAME);

    if (!report_hist_from_read_stats(length_stats, &length)
        || !report_hist_from_read_stats(quality_stats, &quality)
        || length.total_count <= 0.0 || quality.total_count <= 0.0) {
        FC_LOG_ERROR("Cannot write HTML report: no in-memory length/quality histograms for '%s'.", output_dir);
        goto cleanup;
    }

    char *data_json = build_report_data_json(&length, &quality, length_qscore, barcode_index);
    char *html = render_report_html(data_json, sample, sample_type, barcode, run_config, run_config_count);
    if (html && write_text_file(report_path, html) == 0) {
        FC_LOG_INFO("Wrote HTML report '%s'.", report_path);
        rc = 0;
    }
    free(data_json);
    free(html);

cleanup:
    report_hist_destroy(&length);
    report_hist_destroy(&quality);
    free(report_path);
    return rc;
}

int fc_write_html_reports(const fc_report_opts *opts)
{
    if (!opts || !opts->output_root) {
        return -1;
    }

    const char *sample_type = "fastcat";
    if (!opts->use_barcodes) {
        if (!opts->engines || opts->engine_count == 0) {
            FC_LOG_ERROR("Cannot write HTML report: no stats engine available.");
            return -1;
        }
        const read_stats *length = NULL;
        const read_stats *quality = NULL;
        select_report_stats(opts->engines[0], opts->unmapped_override, &length, &quality);
        return write_report_for_output_dir(
            opts->output_root,
            opts->sample,
            sample_type,
            "unclassified",
            length,
            quality,
            opts->run_config,
            opts->run_config_count,
            opts->length_qscore,
            0);
    }

    int status = 0;
    size_t barcode_count = opts->barcode_count ? opts->barcode_count : 1;
    if (!opts->engines || opts->engine_count == 0) {
        FC_LOG_ERROR("Cannot write HTML reports: no stats engines available.");
        return -1;
    }
    if (barcode_count > opts->engine_count) {
        barcode_count = opts->engine_count;
    }
    for (size_t barcode = 0; barcode < barcode_count; ++barcode) {
        char barcode_buf[32] = { 0 };
        const char *barcode_name = fc_sink_barcode_stem(barcode, barcode_buf, sizeof(barcode_buf));
        if (!barcode_name) {
            status = -1;
            continue;
        }
        const read_stats *length = NULL;
        const read_stats *quality = NULL;
        select_report_stats(opts->engines[barcode], opts->unmapped_override, &length, &quality);
        bool active = read_stats_total_count(length) > 0.0 && read_stats_total_count(quality) > 0.0;
        if (!active) {
            continue;
        }
        char *barcode_dir = fc_join_path(opts->output_root, barcode_name);
        if (write_report_for_output_dir(
                barcode_dir,
                opts->sample ? opts->sample : barcode_name,
                sample_type,
                barcode_name,
                length,
                quality,
                opts->run_config,
                opts->run_config_count,
                opts->length_qscore,
                barcode)
            != 0) {
            status = -1;
        }
        free(barcode_dir);
    }
    return status;
}
