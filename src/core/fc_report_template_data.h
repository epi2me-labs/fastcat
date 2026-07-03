#ifndef FC_REPORT_TEMPLATE_HTML_DATA_H
#define FC_REPORT_TEMPLATE_HTML_DATA_H

#include <stddef.h>

typedef struct {
    const char *name;
    const char *model_id;
    const char *attr;
} fc_report_template_binding;

extern const unsigned char fc_report_template_html_gz[];
extern const size_t fc_report_template_html_gz_len;
extern const size_t fc_report_template_html_len;

extern const unsigned char fc_report_template_html_docs_json_gz[];
extern const size_t fc_report_template_html_docs_json_gz_len;
extern const size_t fc_report_template_html_docs_json_len;

extern const fc_report_template_binding fc_report_template_html_bindings[];
extern const size_t fc_report_template_html_binding_count;

#endif
