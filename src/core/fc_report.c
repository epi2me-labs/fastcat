#include "core/fc_report.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "htslib/kstring.h"

#include "core/fc_common.h"
#include "core/fc_json.h"
#include "core/fc_log.h"
#include "core/fc_report_template_data.h"
#include "sink/common.h"
#include "vendor/parson/parson.h"
#include "version.h"

#define FC_REPORT_LENGTH_BINS 100
#define FC_REPORT_QUALITY_MIN 4.0
#define FC_REPORT_QUALITY_MAX 30.0
#define FC_REPORT_QUALITY_WIDTH 0.2
#define FC_REPORT_QUALITY_BINS 130
#define FC_REPORT_YIELD_THINNING 100
#define FC_REPORT_LQ2D_COLOR_QUANTILE 0.995
#define FC_REPORT_TIME_BOX_MAX_BINS 24
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

typedef struct {
    fc_json_value_t *value;
    fc_json_array_t *entries;
} report_source;

typedef struct {
    fc_json_value_t *value;
    fc_json_array_t *array;
} report_column;

// string helpers
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

//
// read stats conversion
//

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

//
// html template helpers
//

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

//
// bokeh document mutation helpers
//

static fc_json_value_t *bokeh_source_data_value(fc_json_array_t **entries)
{
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);
    fc_json_object_set_string(object, "type", "map");
    *entries = fc_json_object_set_array(object, "entries");
    return value;
}

static report_source bokeh_source_create(void)
{
    report_source source = { 0 };
    source.value = bokeh_source_data_value(&source.entries);
    return source;
}

static report_column bokeh_source_add_array_column(report_source *source, const char *name)
{
    report_column column = { 0 };
    column.value = fc_json_array_value();
    column.array = fc_json_get_array(column.value);
    fc_json_entries_append_value(source->entries, name, column.value);
    return column;
}

static JSON_Object *bokeh_find_model(JSON_Value *value, const char *model_id)
{
    if (!value || !model_id) {
        return NULL;
    }

    JSON_Value_Type type = json_value_get_type(value);
    if (type == JSONObject) {
        JSON_Object *object = json_value_get_object(value);
        const char *id = json_object_get_string(object, "id");
        if (id && strcmp(id, model_id) == 0) {
            return object;
        }

        size_t count = json_object_get_count(object);
        for (size_t i = 0; i < count; i++) {
            JSON_Object *hit = bokeh_find_model(json_object_get_value_at(object, i), model_id);
            if (hit) {
                return hit;
            }
        }
        return NULL;
    }

    if (type == JSONArray) {
        JSON_Array *array = json_value_get_array(value);
        size_t count = json_array_get_count(array);
        for (size_t i = 0; i < count; i++) {
            JSON_Object *hit = bokeh_find_model(json_array_get_value(array, i), model_id);
            if (hit) {
                return hit;
            }
        }
    }

    return NULL;
}

static void bokeh_set_value(JSON_Value *docs, const char *name, fc_json_value_t *value)
{
    const fc_report_template_binding *binding = NULL;
    for (size_t i = 0; i < fc_report_template_html_binding_count; i++) {
        if (strcmp(fc_report_template_html_bindings[i].name, name) == 0) {
            binding = &fc_report_template_html_bindings[i];
            break;
        }
    }
    if (!binding) {
        fc_json_value_free(value);
        FC_LOG_FATAL("Report template binding '%s' is missing.", name);
    }

    JSON_Object *model = bokeh_find_model(docs, binding->model_id);
    if (!model) {
        fc_json_value_free(value);
        FC_LOG_FATAL("Bokeh model '%s' for binding '%s' is missing.", binding->model_id, name);
    }

    JSON_Object *attrs = json_object_get_object(model, "attributes");
    if (!attrs) {
        fc_json_value_free(value);
        FC_LOG_FATAL("Bokeh model '%s' has no attributes object.", binding->model_id);
    }

    if (json_object_set_value(attrs, binding->attr, value) != JSONSuccess) {
        FC_LOG_FATAL("Failed to update Bokeh binding '%s'.", name);
    }
}

static void bokeh_set_string(JSON_Value *docs, const char *name, const char *text)
{
    bokeh_set_value(docs, name, fc_json_string_value(text));
}

static void bokeh_set_number(JSON_Value *docs, const char *name, double value)
{
    bokeh_set_value(docs, name, fc_json_number_value(value));
}

static void bokeh_set_range(JSON_Value *docs, const char *prefix, double start, double end)
{
    char name[96];
    double bounds[2] = { start, end };
    FC_SNPRINTF_CHECK(name, sizeof(name), "%s.start", prefix);
    bokeh_set_number(docs, name, start);
    FC_SNPRINTF_CHECK(name, sizeof(name), "%s.end", prefix);
    bokeh_set_number(docs, name, end);
    FC_SNPRINTF_CHECK(name, sizeof(name), "%s.bounds", prefix);
    bokeh_set_value(docs, name, fc_json_number_array_value(bounds, 2));
}

//
// embedded report template helpers
//

static char *inflate_report_asset(
    const unsigned char *compressed,
    size_t compressed_len,
    size_t raw_len,
    const char *label)
{
    if (raw_len > SIZE_MAX - 1 || raw_len > UINT_MAX || compressed_len > UINT_MAX) {
        FC_LOG_ERROR("Embedded report %s is too large to inflate.", label);
        return NULL;
    }

    char *out = fc_alloc(raw_len + 1, sizeof(*out), label);
    z_stream stream = { 0 };
    stream.next_in = (Bytef *)compressed;
    stream.avail_in = (uInt)compressed_len;
    stream.next_out = (Bytef *)out;
    stream.avail_out = (uInt)raw_len;

    int zret = inflateInit2(&stream, 16 + MAX_WBITS);
    if (zret != Z_OK) {
        FC_LOG_ERROR("Failed to initialise report %s inflater.", label);
        free(out);
        return NULL;
    }
    zret = inflate(&stream, Z_FINISH);
    int zend = inflateEnd(&stream);

    if (zret != Z_STREAM_END || zend != Z_OK || stream.total_out != raw_len) {
        FC_LOG_ERROR("Failed to inflate embedded report %s.", label);
        free(out);
        return NULL;
    }

    out[raw_len] = '\0';
    return out;
}

//
// histogram summary helpers
//

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

//
// quality plot
//

static fc_json_value_t *build_quality_source_data_value(
    const report_hist *quality,
    char *subtitle,
    size_t subtitle_size,
    double *y_end)
{
    double top[FC_REPORT_QUALITY_BINS] = { 0 };
    double left[FC_REPORT_QUALITY_BINS] = { 0 };
    double right[FC_REPORT_QUALITY_BINS] = { 0 };
    *y_end = 1.0;

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
            if (top[bin] > *y_end) {
                *y_end = top[bin];
            }
        }
    }
    *y_end *= 1.1;

    FC_SNPRINTF_CHECK(subtitle, subtitle_size, "Mean: %.1f. Median: %.1f",
        hist_mean_midpoint(quality), hist_median_start(quality));

    report_source source = bokeh_source_create();
    fc_json_entries_append_value(
        source.entries, "top", fc_json_number_array_value(top, FC_REPORT_QUALITY_BINS));
    fc_json_entries_append_value(
        source.entries, "left", fc_json_number_array_value(left, FC_REPORT_QUALITY_BINS));
    fc_json_entries_append_value(
        source.entries, "right", fc_json_number_array_value(right, FC_REPORT_QUALITY_BINS));
    return source.value;
}

static void update_quality_plot(JSON_Value *docs, const report_hist *quality)
{
    char subtitle[128] = { 0 };
    double y_end = 1.0;
    bokeh_set_value(
        docs,
        "quality.source.data",
        build_quality_source_data_value(quality, subtitle, sizeof(subtitle), &y_end));
    bokeh_set_string(docs, "quality.subtitle.text", subtitle);
    bokeh_set_range(docs, "quality.x_range", FC_REPORT_QUALITY_MIN, FC_REPORT_QUALITY_MAX);
    bokeh_set_range(docs, "quality.y_range", 0.0, y_end);
}

//
// length plot
//

static fc_json_value_t *build_length_source_data_value(
    const report_hist *length,
    char *subtitle,
    size_t subtitle_size,
    double *x_start,
    double *x_end,
    double *y_end)
{
    double top[FC_REPORT_LENGTH_BINS] = { 0 };
    double left[FC_REPORT_LENGTH_BINS] = { 0 };
    double right[FC_REPORT_LENGTH_BINS] = { 0 };

    double min_kb = (length->n > 0) ? length->bins[0].start / 1000.0 : 0.0;
    double max_kb = (length->n > 0) ? length->bins[length->n - 1].start / 1000.0 : 1.0;
    if (max_kb <= min_kb) {
        max_kb = min_kb + 1.0;
    }
    *x_start = min_kb;
    *x_end = max_kb;
    *y_end = 1.0;

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
            if (top[bin] > *y_end) {
                *y_end = top[bin];
            }
        }
    }
    *y_end *= 1.1;

    uint64_t min_len = (length->n > 0) ? (uint64_t)llround(length->bins[0].start) : 0;
    uint64_t max_len = (length->n > 0) ? (uint64_t)llround(length->bins[length->n - 1].start) : 0;
    FC_SNPRINTF_CHECK(subtitle, subtitle_size,
        "Mean: %.1f. Median: %.1f. Min: %" PRIu64 ". Max: %" PRIu64,
        hist_mean_start(length), hist_median_start(length), min_len, max_len);

    report_source source = bokeh_source_create();
    fc_json_entries_append_value(
        source.entries, "top", fc_json_number_array_value(top, FC_REPORT_LENGTH_BINS));
    fc_json_entries_append_value(
        source.entries, "left", fc_json_number_array_value(left, FC_REPORT_LENGTH_BINS));
    fc_json_entries_append_value(
        source.entries, "right", fc_json_number_array_value(right, FC_REPORT_LENGTH_BINS));
    return source.value;
}

static void update_length_plot(JSON_Value *docs, const report_hist *length)
{
    char subtitle[192] = { 0 };
    double x_start = 0.0;
    double x_end = 1.0;
    double y_end = 1.0;
    bokeh_set_value(
        docs,
        "length.source.data",
        build_length_source_data_value(length, subtitle, sizeof(subtitle), &x_start, &x_end, &y_end));
    bokeh_set_string(docs, "length.subtitle.text", subtitle);
    bokeh_set_range(docs, "length.x_range", x_start, x_end);
    bokeh_set_range(docs, "length.y_range", 0.0, y_end);
}

//
// yield plot
//

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

static fc_json_value_t *build_yield_source_data_value(
    const report_hist *length,
    char *subtitle,
    size_t subtitle_size,
    double *x_end,
    double *y_end)
{
    double *x = NULL;
    double *y = NULL;
    size_t out_n = 0;
    size_t out_cap = 0;
    double *cumsum = NULL;
    double total_yield = 0.0;
    *x_end = 1.0;
    *y_end = 1.0;

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
        *x_end = length->bins[length->n - 1].start / 1000.0;
        if (!isfinite(*x_end) || *x_end <= 0.0) {
            *x_end = 1.0;
        }
        *y_end = total_yield / 1e9;
        *y_end = (isfinite(*y_end) && *y_end > 0.0) ? *y_end * 1.1 : 1.0;
    }

    FC_SNPRINTF_CHECK(subtitle, subtitle_size, "Total yield: %.3g Gb.  N50: %.3gkb",
        total_yield / 1e9, compute_n50_kb(length, total_yield));

    report_source source = bokeh_source_create();
    fc_json_entries_append_value(source.entries, "x", fc_json_number_array_value(x, out_n));
    fc_json_entries_append_value(source.entries, "y", fc_json_number_array_value(y, out_n));

    free(cumsum);
    free(x);
    free(y);
    return source.value;
}

static void update_yield_plot(JSON_Value *docs, const report_hist *length)
{
    char subtitle[160] = { 0 };
    double x_end = 1.0;
    double y_end = 1.0;
    bokeh_set_value(
        docs,
        "yield.source.data",
        build_yield_source_data_value(length, subtitle, sizeof(subtitle), &x_end, &y_end));
    bokeh_set_string(docs, "yield.subtitle.text", subtitle);
    bokeh_set_range(docs, "yield.x_range", 0.0, x_end);
    bokeh_set_range(docs, "yield.y_range", 0.0, y_end);
}

//
// length/qscore heatmap
//

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

static fc_json_value_t *build_lq2d_empty_source_data_value(void)
{
    report_source source = bokeh_source_create();
    bokeh_source_add_array_column(&source, "image");
    bokeh_source_add_array_column(&source, "x");
    bokeh_source_add_array_column(&source, "y");
    bokeh_source_add_array_column(&source, "dw");
    bokeh_source_add_array_column(&source, "dh");
    return source.value;
}

static char *build_lq2d_image_base64(const fc_length_qscore_2d_matrix *matrix)
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

    kstring_t b64 = { 0, 0, NULL };
    append_base64(&b64, bytes, byte_count);
    free(bytes);
    return report_steal(&b64);
}

static fc_json_value_t *build_lq2d_image_ndarray_value(const fc_length_qscore_2d_matrix *matrix)
{
    char *b64 = build_lq2d_image_base64(matrix);

    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);
    fc_json_object_set_string(object, "type", "ndarray");

    fc_json_object_t *array = fc_json_object_set_object(object, "array");
    fc_json_object_set_string(array, "type", "bytes");
    fc_json_object_set_string(array, "data", b64);

    fc_json_array_t *shape = fc_json_object_set_array(object, "shape");
    fc_json_array_append_number(shape, (double)matrix->length_bins);
    fc_json_array_append_number(shape, (double)matrix->qscore_bins);

    fc_json_object_set_string(object, "dtype", "float64");
    fc_json_object_set_string(object, "order", "little");

    free(b64);
    return value;
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

static void update_length_qscore_plot(
    JSON_Value *docs,
    const fc_length_qscore_2d *model,
    size_t barcode)
{
    if (!model || !fc_length_qscore_2d_seen(model, barcode)) {
        bokeh_set_value(
            docs,
            "length_qscore.source.data",
            build_lq2d_empty_source_data_value());
        bokeh_set_string(
            docs,
            "length_qscore.subtitle.text",
            "No length/qscore data available");
        bokeh_set_number(docs, "length_qscore.color_mapper.high", 1.0);
        bokeh_set_range(docs, "length_qscore.x_range", 0.0, 1.0);
        bokeh_set_range(docs, "length_qscore.y_range", 0.0, 1.0);
        return;
    }

    fc_length_qscore_2d_matrix matrix = { 0 };
    if (fc_length_qscore_2d_materialize(model, barcode, &matrix) != 0 || !matrix.counts) {
        fc_length_qscore_2d_matrix_destroy(&matrix);
        update_length_qscore_plot(docs, NULL, barcode);
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
    double x_start = 0.0;
    double y_start = 0.0;
    double x_end = matrix.qscore_width * (double)matrix.qscore_bins;
    double y_end = matrix.length_upper / 1000.0;

    report_source source = bokeh_source_create();
    report_column image = bokeh_source_add_array_column(&source, "image");
    report_column x = bokeh_source_add_array_column(&source, "x");
    report_column y = bokeh_source_add_array_column(&source, "y");
    report_column dw = bokeh_source_add_array_column(&source, "dw");
    report_column dh = bokeh_source_add_array_column(&source, "dh");
    fc_json_array_append_value(image.array, build_lq2d_image_ndarray_value(&matrix));
    fc_json_array_append_number(x.array, x_start);
    fc_json_array_append_number(y.array, y_start);
    fc_json_array_append_number(dw.array, x_end - x_start);
    fc_json_array_append_number(dh.array, y_end - y_start);

    bokeh_set_value(docs, "length_qscore.source.data", source.value);
    bokeh_set_string(docs, "length_qscore.subtitle.text", subtitle);
    bokeh_set_number(
        docs, "length_qscore.color_mapper.high", color_high > 0.0 ? color_high : 1.0);
    bokeh_set_range(docs, "length_qscore.x_range", x_start, x_end);
    bokeh_set_range(docs, "length_qscore.y_range", y_start, y_end);

    fc_length_qscore_2d_matrix_destroy(&matrix);
}

//
// time plots
//

static void format_report_hours(double hours, char *buf, size_t size)
{
    double rounded = round(hours);
    if (fabs(hours - rounded) < 1e-9) {
        FC_SNPRINTF_CHECK(buf, size, "%.0f", rounded);
        return;
    }

    FC_SNPRINTF_CHECK(buf, size, "%.2f", hours);
    char *dot = strchr(buf, '.');
    if (!dot) {
        return;
    }

    char *end = buf + strlen(buf) - 1;
    while (end > dot && *end == '0') {
        *end = '\0';
        end--;
    }
    if (end == dot) {
        *end = '\0';
    }
}

static void format_time_bin_label(const fc_time_bin_row *row, char *buf, size_t size)
{
    char start[32] = { 0 };
    char end[32] = { 0 };
    format_report_hours(row->time_bin_start_s / 3600.0, start, sizeof(start));
    format_report_hours(row->time_bin_end_s / 3600.0, end, sizeof(end));
    FC_SNPRINTF_CHECK(buf, size, "%s-%s", start, end);
}

static bool time_quantiles_valid(const fc_time_bin_quantiles *q)
{
    return q && isfinite(q->q05) && isfinite(q->q25) && isfinite(q->q50)
        && isfinite(q->q75) && isfinite(q->q95);
}

static bool time_box_row_valid(const fc_time_bin_row *row, bool length)
{
    return row && row->read_count > 0
        && time_quantiles_valid(length ? &row->length_quantiles : &row->qscore_quantiles);
}

static void time_box_y_range(
    const fc_time_bin_table *table,
    bool length,
    double *start,
    double *end)
{
    double y_min = INFINITY;
    double y_max = -INFINITY;
    for (size_t i = 0; table && i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        if (!time_box_row_valid(row, length)) {
            continue;
        }
        const fc_time_bin_quantiles *q = length ? &row->length_quantiles : &row->qscore_quantiles;
        double scale = length ? 1000.0 : 1.0;
        double lower = q->q05 / scale;
        double upper = q->q95 / scale;
        if (lower < y_min) {
            y_min = lower;
        }
        if (upper > y_max) {
            y_max = upper;
        }
    }

    if (!isfinite(y_min) || !isfinite(y_max)) {
        *start = 0.0;
        *end = 1.0;
        return;
    }

    double padding = (y_max - y_min) * 0.1;
    if (padding == 0.0) {
        padding = fabs(y_max) * 0.1;
    }
    if (padding == 0.0) {
        padding = 1.0;
    }

    *start = y_min - padding;
    if (y_min >= 0.0 && *start < 0.0) {
        *start = 0.0;
    }
    *end = y_max + padding;
    if (*end <= *start) {
        *end = *start + 1.0;
    }
}

static void update_time_yield_plot(JSON_Value *docs, const fc_time_bin_table *table)
{
    report_source source = bokeh_source_create();
    report_column x = bokeh_source_add_array_column(&source, "x");
    report_column y = bokeh_source_add_array_column(&source, "y");
    double y_end = 0.0;

    for (size_t i = 0; table && i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        double mid_h = (row->time_bin_start_s + row->time_bin_end_s) / 7200.0;
        y_end = (double)row->cumulative_base_count / 1e9;
        fc_json_array_append_number(x.array, mid_h);
        fc_json_array_append_number(y.array, y_end);
    }

    bokeh_set_value(docs, "time_yield.source.data", source.value);

    double x_end = (table && table->n > 0) ? table->rows[table->n - 1].time_bin_end_s / 3600.0 : 1.0;
    if (!isfinite(x_end) || x_end <= 0.0) {
        x_end = 1.0;
    }
    if (!isfinite(y_end) || y_end <= 0.0) {
        y_end = 1.0;
    } else {
        y_end *= 1.1;
    }
    bokeh_set_range(docs, "time_yield.x_range", 0.0, x_end);
    bokeh_set_range(docs, "time_yield.y_range", 0.0, y_end);
}

static void update_time_reads_plot(JSON_Value *docs, const fc_time_bin_table *table)
{
    report_source source = bokeh_source_create();
    report_column x = bokeh_source_add_array_column(&source, "x");
    report_column y = bokeh_source_add_array_column(&source, "y");
    double y_end = 0.0;

    for (size_t i = 0; table && i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        double mid_h = (row->time_bin_start_s + row->time_bin_end_s) / 7200.0;
        y_end = (double)row->cumulative_read_count;
        fc_json_array_append_number(x.array, mid_h);
        fc_json_array_append_number(y.array, y_end);
    }

    bokeh_set_value(docs, "time_reads.source.data", source.value);

    double x_end = (table && table->n > 0) ? table->rows[table->n - 1].time_bin_end_s / 3600.0 : 1.0;
    if (!isfinite(x_end) || x_end <= 0.0) {
        x_end = 1.0;
    }
    if (!isfinite(y_end) || y_end <= 0.0) {
        y_end = 1.0;
    } else {
        y_end *= 1.1;
    }
    bokeh_set_range(docs, "time_reads.x_range", 0.0, x_end);
    bokeh_set_range(docs, "time_reads.y_range", 0.0, y_end);
}

static void update_time_length_box_plot(JSON_Value *docs, const fc_time_bin_table *table)
{
    report_source source = bokeh_source_create();
    report_column index = bokeh_source_add_array_column(&source, "index");
    report_column labels = bokeh_source_add_array_column(&source, "time bin (hours)");
    report_column read_count = bokeh_source_add_array_column(&source, "read_count");
    report_column q05 = bokeh_source_add_array_column(&source, "q05");
    report_column q25 = bokeh_source_add_array_column(&source, "q25");
    report_column q50 = bokeh_source_add_array_column(&source, "read length (kb)");
    report_column q75 = bokeh_source_add_array_column(&source, "q75");
    report_column q95 = bokeh_source_add_array_column(&source, "q95");
    report_column mean = bokeh_source_add_array_column(&source, "mean");
    report_column cap_left = bokeh_source_add_array_column(&source, "cap_left");
    report_column cap_right = bokeh_source_add_array_column(&source, "cap_right");
    fc_json_value_t *factors = fc_json_array_value();
    fc_json_array_t *factors_array = fc_json_get_array(factors);

    size_t out_i = 0;
    for (size_t i = 0; table && i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        if (!time_box_row_valid(row, true)) {
            continue;
        }

        const fc_time_bin_quantiles *q = &row->length_quantiles;
        double scaled_q05 = q->q05 / 1000.0;
        double scaled_q25 = q->q25 / 1000.0;
        double scaled_q50 = q->q50 / 1000.0;
        double scaled_q75 = q->q75 / 1000.0;
        double scaled_q95 = q->q95 / 1000.0;
        double scaled_mean = (isfinite(row->length_mean) ? row->length_mean : q->q50) / 1000.0;
        char label[64] = { 0 };
        format_time_bin_label(row, label, sizeof(label));

        fc_json_array_append_number(index.array, (double)out_i);
        fc_json_array_append_string(labels.array, label);
        fc_json_array_append_string(factors_array, label);
        fc_json_array_append_number(read_count.array, (double)row->read_count);
        fc_json_array_append_number(q05.array, scaled_q05);
        fc_json_array_append_number(q25.array, scaled_q25);
        fc_json_array_append_number(q50.array, scaled_q50);
        fc_json_array_append_number(q75.array, scaled_q75);
        fc_json_array_append_number(q95.array, scaled_q95);
        fc_json_array_append_number(mean.array, scaled_mean);

        fc_json_array_t *left = fc_json_array_append_array(cap_left.array);
        fc_json_array_append_string(left, label);
        fc_json_array_append_number(left, -0.15);
        fc_json_array_t *right = fc_json_array_append_array(cap_right.array);
        fc_json_array_append_string(right, label);
        fc_json_array_append_number(right, 0.15);

        out_i++;
    }

    double y_start = 0.0;
    double y_end = 1.0;
    time_box_y_range(table, true, &y_start, &y_end);
    bokeh_set_value(docs, "time_length.source.data", source.value);
    bokeh_set_value(docs, "time_length.x_range.factors", factors);
    bokeh_set_range(docs, "time_length.y_range", y_start, y_end);
}

static void update_time_qscore_box_plot(JSON_Value *docs, const fc_time_bin_table *table)
{
    report_source source = bokeh_source_create();
    report_column index = bokeh_source_add_array_column(&source, "index");
    report_column labels = bokeh_source_add_array_column(&source, "time bin (hours)");
    report_column read_count = bokeh_source_add_array_column(&source, "read_count");
    report_column q05 = bokeh_source_add_array_column(&source, "q05");
    report_column q25 = bokeh_source_add_array_column(&source, "q25");
    report_column q50 = bokeh_source_add_array_column(&source, "qscore");
    report_column q75 = bokeh_source_add_array_column(&source, "q75");
    report_column q95 = bokeh_source_add_array_column(&source, "q95");
    report_column mean = bokeh_source_add_array_column(&source, "mean");
    report_column cap_left = bokeh_source_add_array_column(&source, "cap_left");
    report_column cap_right = bokeh_source_add_array_column(&source, "cap_right");
    fc_json_value_t *factors = fc_json_array_value();
    fc_json_array_t *factors_array = fc_json_get_array(factors);

    size_t out_i = 0;
    for (size_t i = 0; table && i < table->n; i++) {
        const fc_time_bin_row *row = &table->rows[i];
        if (!time_box_row_valid(row, false)) {
            continue;
        }

        const fc_time_bin_quantiles *q = &row->qscore_quantiles;
        double scaled_mean = isfinite(row->qscore_mean) ? row->qscore_mean : q->q50;
        char label[64] = { 0 };
        format_time_bin_label(row, label, sizeof(label));

        fc_json_array_append_number(index.array, (double)out_i);
        fc_json_array_append_string(labels.array, label);
        fc_json_array_append_string(factors_array, label);
        fc_json_array_append_number(read_count.array, (double)row->read_count);
        fc_json_array_append_number(q05.array, q->q05);
        fc_json_array_append_number(q25.array, q->q25);
        fc_json_array_append_number(q50.array, q->q50);
        fc_json_array_append_number(q75.array, q->q75);
        fc_json_array_append_number(q95.array, q->q95);
        fc_json_array_append_number(mean.array, scaled_mean);

        fc_json_array_t *left = fc_json_array_append_array(cap_left.array);
        fc_json_array_append_string(left, label);
        fc_json_array_append_number(left, -0.15);
        fc_json_array_t *right = fc_json_array_append_array(cap_right.array);
        fc_json_array_append_string(right, label);
        fc_json_array_append_number(right, 0.15);

        out_i++;
    }

    double y_start = 0.0;
    double y_end = 1.0;
    time_box_y_range(table, false, &y_start, &y_end);
    bokeh_set_value(docs, "time_qscore.source.data", source.value);
    bokeh_set_value(docs, "time_qscore.x_range.factors", factors);
    bokeh_set_range(docs, "time_qscore.y_range", y_start, y_end);
}

//
// report assembly
//

static JSON_Value *build_report_bokeh_docs(
    const report_hist *length,
    const report_hist *quality,
    const fc_length_qscore_2d *length_qscore,
    const fc_time_bin *time_bin,
    size_t barcode)
{
    char *docs_json = inflate_report_asset(
        fc_report_template_html_docs_json_gz,
        fc_report_template_html_docs_json_gz_len,
        fc_report_template_html_docs_json_len,
        "Bokeh document JSON");
    if (!docs_json) {
        return NULL;
    }

    JSON_Value *docs = fc_json_parse_string(docs_json);
    free(docs_json);

    fc_time_bin_table time_table = { 0 };
    if (time_bin && fc_time_bin_materialize(time_bin, barcode, &time_table) != 0) {
        fc_time_bin_table_destroy(&time_table);
        memset(&time_table, 0, sizeof(time_table));
    }

    fc_time_bin_table time_boxes = { 0 };
    if (time_bin
        && fc_time_bin_materialize_merged(time_bin, barcode, FC_REPORT_TIME_BOX_MAX_BINS, &time_boxes)
            != 0) {
        fc_time_bin_table_destroy(&time_boxes);
        memset(&time_boxes, 0, sizeof(time_boxes));
    }

    update_quality_plot(docs, quality);
    update_length_plot(docs, length);
    update_yield_plot(docs, length);
    update_length_qscore_plot(docs, length_qscore, barcode);
    update_time_yield_plot(docs, &time_table);
    update_time_reads_plot(docs, &time_table);
    update_time_length_box_plot(docs, &time_boxes);
    update_time_qscore_box_plot(docs, &time_boxes);

    fc_time_bin_table_destroy(&time_boxes);
    fc_time_bin_table_destroy(&time_table);
    return docs;
}

//
// report rendering and output
//

static char *render_report_html(
    const char *docs_json,
    const char *sample,
    const char *sample_type,
    const char *barcode,
    const fc_report_config_entry *run_config,
    size_t run_config_count)
{
    char *html = inflate_report_asset(
        fc_report_template_html_gz,
        fc_report_template_html_gz_len,
        fc_report_template_html_len,
        "HTML template");
    if (!html) {
        return NULL;
    }

    char *sample_html = html_escape(sample ? sample : "fastcat");
    char *type_html = html_escape(sample_type ? sample_type : "fastcat");
    char *barcode_html = html_escape(barcode ? barcode : "unclassified");
    char *version_html = html_escape(argp_program_version ? argp_program_version : "unknown");
    char *run_config_rows_html = build_run_config_rows_html(run_config, run_config_count);

    replace_all(&html, "__FASTCAT_BOKEH_DOCS_JSON__", docs_json ? docs_json : "{}");
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
    const fc_time_bin *time_bin,
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

    JSON_Value *docs = build_report_bokeh_docs(&length, &quality, length_qscore, time_bin, barcode_index);
    char *docs_json = docs ? fc_json_serialize_to_string(docs) : NULL;
    char *html = render_report_html(docs_json, sample, sample_type, barcode, run_config, run_config_count);
    if (html && write_text_file(report_path, html) == 0) {
        FC_LOG_INFO("Wrote HTML report '%s'.", report_path);
        rc = 0;
    }
    if (docs) {
        fc_json_value_free(docs);
    }
    free(docs_json);
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
            opts->time_bin,
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
                opts->time_bin,
                barcode)
            != 0) {
            status = -1;
        }
        free(barcode_dir);
    }
    return status;
}
