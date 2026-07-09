#ifndef FASTCAT_FC_JSON_H
#define FASTCAT_FC_JSON_H

#include <stddef.h>

typedef struct json_value_t fc_json_value_t;
typedef struct json_object_t fc_json_object_t;
typedef struct json_array_t fc_json_array_t;

fc_json_value_t *fc_json_object_value(void);
fc_json_value_t *fc_json_array_value(void);
fc_json_value_t *fc_json_string_value(const char *value);
fc_json_value_t *fc_json_number_value(double value);
fc_json_value_t *fc_json_bool_value(int value);
fc_json_value_t *fc_json_null_value(void);
fc_json_value_t *fc_json_parse_string(const char *json);
fc_json_object_t *fc_json_get_object(fc_json_value_t *value);
fc_json_array_t *fc_json_get_array(fc_json_value_t *value);
void fc_json_value_free(fc_json_value_t *value);
char *fc_json_serialize_to_string(const fc_json_value_t *value);

void fc_json_object_set_value(fc_json_object_t *object, const char *key, fc_json_value_t *value);
void fc_json_object_set_string(fc_json_object_t *object, const char *key, const char *value);
void fc_json_object_set_number(fc_json_object_t *object, const char *key, double value);
void fc_json_object_set_null(fc_json_object_t *object, const char *key);
fc_json_object_t *fc_json_object_set_object(fc_json_object_t *object, const char *key);
fc_json_array_t *fc_json_object_set_array(fc_json_object_t *object, const char *key);

void fc_json_array_append_value(fc_json_array_t *array, fc_json_value_t *value);
void fc_json_array_append_string(fc_json_array_t *array, const char *value);
void fc_json_array_append_number(fc_json_array_t *array, double value);
fc_json_array_t *fc_json_array_append_array(fc_json_array_t *array);

fc_json_value_t *fc_json_number_array_value(const double *values, size_t n);
void fc_json_entries_append_value(fc_json_array_t *entries, const char *key, fc_json_value_t *value);

#endif
