#include "core/fc_report_arg_table.h"

#include <inttypes.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/kstring.h"

#include "core/fc_common.h"
#include "core/fc_log.h"

typedef struct {
    fc_report_config_entry *entries;
    size_t n_entries;
    size_t cap;
} report_config_builder;

static char *join_string_values(const char *const *values, size_t n)
{
    if (!values || n == 0) {
        return fc_strdup("(none)", "report configuration value");
    }

    kstring_t value = { 0, 0, NULL };
    for (size_t i = 0; i < n; ++i) {
        if (i > 0 && kputs(", ", &value) < 0) {
            free(value.s);
            FC_LOG_FATAL("Failed to format report configuration value.");
        }
        if (kputs(values[i] ? values[i] : "(not set)", &value) < 0) {
            free(value.s);
            FC_LOG_FATAL("Failed to format report configuration value.");
        }
    }
    return value.s;
}

static char *join_u32_values(const uint32_t *values, size_t n)
{
    if (!values || n == 0) {
        return fc_strdup("(none)", "report configuration value");
    }

    kstring_t value = { 0, 0, NULL };
    for (size_t i = 0; i < n; ++i) {
        if (i > 0 && kputs(", ", &value) < 0) {
            free(value.s);
            FC_LOG_FATAL("Failed to format report configuration value.");
        }
        if (ksprintf(&value, "%" PRIu32, values[i]) < 0) {
            free(value.s);
            FC_LOG_FATAL("Failed to format report configuration value.");
        }
    }
    return value.s;
}

static void config_add_str(report_config_builder *builder, const char *key, const char *value)
{
    if (builder->n_entries == builder->cap) {
        size_t new_cap = builder->cap ? builder->cap * 2 : 32;
        builder->entries = fc_realloc(
            builder->entries,
            new_cap * sizeof(*builder->entries),
            "report configuration entries");
        builder->cap = new_cap;
    }
    builder->entries[builder->n_entries++] = (fc_report_config_entry) {
        .key = key,
        .value = fc_strdup(value ? value : "(not set)", "report configuration value")
    };
}

static void config_add_bool(report_config_builder *builder, const char *key, bool value)
{
    config_add_str(builder, key, value ? "true" : "false");
}

static void config_add_strf(report_config_builder *builder, const char *key, const char *fmt, ...)
{
    kstring_t value = { 0, 0, NULL };
    va_list ap;
    va_start(ap, fmt);
    int rc = kvsprintf(&value, fmt, ap);
    va_end(ap);
    if (rc < 0) {
        free(value.s);
        FC_LOG_FATAL("Failed to format report configuration value.");
    }
    config_add_str(builder, key, value.s);
    free(value.s);
}

static void config_add_limit(report_config_builder *builder, const char *key, size_t value, bool unlimited)
{
    if (unlimited && value == (size_t)-1) {
        config_add_str(builder, key, "unlimited");
    } else {
        config_add_strf(builder, key, "%zu", value);
    }
}

static const char *config_file_kind(fc_file_kind kind)
{
    switch (kind) {
    case FC_FILE_FASTQ:
        return "fastq";
    case FC_FILE_BAM:
        return "bam";
    default:
        return "unknown";
    }
}

fc_report_arg_table fc_report_arg_table_from_catstats(
    fc_program_kind program,
    const fc_catstats_args *args,
    fc_file_kind input_kind,
    size_t file_count)
{
    report_config_builder builder = { 0 };
    const char *program_name = "fastcat fastq";
    if (program == FC_PROGRAM_BAMSTATS) {
        program_name = "fastcat bam";
    } else if (program == FC_PROGRAM_QC) {
        program_name = "fastcat qc";
    }

    config_add_str(&builder, "subcommand", program_name);
    config_add_str(&builder, "sample", args->sample);
    config_add_str(&builder, "output", args->output_dir);
    config_add_strf(&builder, "input files discovered", "%zu", file_count);
    config_add_bool(&builder, "demultiplex", args->demultiplex);
    config_add_bool(&builder, "per-read output", args->per_read);
    config_add_bool(&builder, "metrics only", args->metrics_only);
    config_add_bool(&builder, "html report", args->html_report);
    config_add_strf(&builder, "threads", "%d", args->threads);
    config_add_bool(&builder, "force error", args->force_error);
    config_add_bool(&builder, "verbose", args->verbose);
    config_add_str(&builder, "log file", args->log_file);

    config_add_limit(&builder, "minimum length", args->min_length, false);
    config_add_limit(&builder, "maximum length", args->max_length, true);
    config_add_strf(&builder, "minimum qscore", "%.3g", args->min_qscore);
    config_add_bool(&builder, "dust filtering", args->dust);
    config_add_strf(&builder, "maximum dust", "%.3g", args->max_dust);
    config_add_limit(&builder, "dust window", args->dust_w, false);
    config_add_limit(&builder, "dust threshold", args->dust_t, false);
    config_add_str(&builder, "read group filter", args->read_group);
    if (args->tag_name[0]) {
        char tag_name[3] = { args->tag_name[0], args->tag_name[1], '\0' };
        config_add_str(&builder, "tag filter name", tag_name);
        config_add_strf(&builder, "tag filter value", "%d", args->tag_value);
    } else {
        config_add_str(&builder, "tag filter name", NULL);
        config_add_str(&builder, "tag filter value", NULL);
    }
    config_add_bool(&builder, "force quality recalculation", args->force_recalc_qual);
    config_add_strf(&builder, "time bin size seconds", "%.10g", args->time_bin_size_s);
    config_add_strf(&builder, "time bin sketch size", "%zu", args->time_bin_k);

    char *input_paths = join_string_values(
        (const char *const *)args->input_paths.paths, args->input_paths.n_paths);
    config_add_str(&builder, "input paths", input_paths);
    free(input_paths);
    config_add_str(&builder, "input format", config_file_kind(input_kind));
    if (program == FC_PROGRAM_FASTCAT || program == FC_PROGRAM_QC) {
        config_add_bool(&builder, "recursive directory search", args->recurse < 0);
        config_add_bool(&builder, "file-of-filenames input", args->fofn);
        config_add_bool(&builder, "BAM output", args->write_bam != 0);
        config_add_limit(&builder, "reads per file", args->reads_per_file, false);
    }
    if (input_kind == FC_FILE_BAM) {
        config_add_str(&builder, "region", args->region);
        config_add_str(&builder, "BED", args->bed);
        config_add_bool(&builder, "coverage", args->coverage);
        char *coverage_beds = join_string_values(args->coverage_beds, args->n_coverage_beds);
        config_add_str(&builder, "coverage BEDs", coverage_beds);
        free(coverage_beds);
        char *coverage_names = join_string_values(args->coverage_names, args->n_coverage_names);
        config_add_str(&builder, "coverage names", coverage_names);
        free(coverage_names);
        char *coverage_thresholds = join_u32_values(args->coverage_thresholds, args->n_coverage_thresholds);
        config_add_str(&builder, "coverage thresholds", coverage_thresholds);
        free(coverage_thresholds);
        char *coverage_segments = join_u32_values(args->segments, args->n_segments);
        config_add_str(&builder, "coverage segments", coverage_segments);
        free(coverage_segments);
        config_add_bool(&builder, "poly-A", args->poly_a);
        config_add_strf(&builder, "poly-A coverage", "%.3g", args->poly_a_cover);
        config_add_strf(&builder, "poly-A qscore", "%.3g", args->poly_a_qual);
        config_add_bool(&builder, "poly-A reverse alignments", args->poly_a_rev);
    }

    return (fc_report_arg_table) {
        .entries = builder.entries,
        .n_entries = builder.n_entries
    };
}

void fc_report_arg_table_destroy(fc_report_arg_table *table)
{
    if (!table) {
        return;
    }
    for (size_t i = 0; i < table->n_entries; ++i) {
        free((char *)table->entries[i].value);
    }
    free(table->entries);
    table->entries = NULL;
    table->n_entries = 0;
}
