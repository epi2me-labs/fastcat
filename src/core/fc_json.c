#include "core/fc_json.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "core/fc_log.h"
#include "vendor/parson/parson.h"

static void *fc_json_malloc(size_t size)
{
    void *ptr = malloc(size ? size : 1);
    if (!ptr) {
        FC_LOG_FATAL("Failed to allocate JSON memory.");
    }
    return ptr;
}

static void fc_json_free(void *ptr)
{
    free(ptr);
}

static void fc_json_init(void)
{
    static bool initialized = false;
    if (!initialized) {
        json_set_allocation_functions(fc_json_malloc, fc_json_free);
        initialized = true;
    }
}

static void fc_json_check(JSON_Status status)
{
    if (status != JSONSuccess) {
        FC_LOG_FATAL("Failed to build JSON.");
    }
}

static double fc_json_number(double value)
{
    return isfinite(value) ? value : 0.0;
}

fc_json_value_t *fc_json_object_value(void)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_object();
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON object.");
    }
    return value;
}

fc_json_value_t *fc_json_array_value(void)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_array();
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON array.");
    }
    return value;
}

fc_json_value_t *fc_json_string_value(const char *text)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_string(text ? text : "");
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON string.");
    }
    return value;
}

fc_json_value_t *fc_json_number_value(double value)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_number(fc_json_number(value));
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON number.");
    }
    return json_value;
}

fc_json_value_t *fc_json_bool_value(int value)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_boolean(value);
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON boolean.");
    }
    return json_value;
}

fc_json_value_t *fc_json_null_value(void)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_null();
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON null.");
    }
    return json_value;
}

fc_json_value_t *fc_json_parse_string(const char *json)
{
    fc_json_init();
    fc_json_value_t *value = json_parse_string(json ? json : "");
    if (!value) {
        FC_LOG_FATAL("Failed to parse JSON.");
    }
    return value;
}

fc_json_object_t *fc_json_get_object(fc_json_value_t *value)
{
    fc_json_object_t *object = json_object(value);
    if (!object) {
        FC_LOG_FATAL("Failed to access JSON object.");
    }
    return object;
}

fc_json_array_t *fc_json_get_array(fc_json_value_t *value)
{
    fc_json_array_t *array = json_array(value);
    if (!array) {
        FC_LOG_FATAL("Failed to access JSON array.");
    }
    return array;
}

void fc_json_value_free(fc_json_value_t *value)
{
    json_value_free(value);
}

char *fc_json_serialize_to_string(const fc_json_value_t *value)
{
    fc_json_init();
    char *serialized = json_serialize_to_string(value);
    if (!serialized) {
        FC_LOG_FATAL("Failed to serialise JSON.");
    }
    return serialized;
}

void fc_json_object_set_value(fc_json_object_t *object, const char *key, fc_json_value_t *value)
{
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON value.");
    }
    fc_json_check(json_object_set_value(object, key, value));
}

void fc_json_object_set_string(fc_json_object_t *object, const char *key, const char *value)
{
    fc_json_check(json_object_set_string(object, key, value ? value : ""));
}

void fc_json_object_set_number(fc_json_object_t *object, const char *key, double value)
{
    fc_json_check(json_object_set_number(object, key, fc_json_number(value)));
}

void fc_json_object_set_null(fc_json_object_t *object, const char *key)
{
    fc_json_check(json_object_set_null(object, key));
}

fc_json_object_t *fc_json_object_set_object(fc_json_object_t *object, const char *key)
{
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_set_value(object, key, value);
    return fc_json_get_object(value);
}

fc_json_array_t *fc_json_object_set_array(fc_json_object_t *object, const char *key)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_object_set_value(object, key, value);
    return fc_json_get_array(value);
}

void fc_json_array_append_value(fc_json_array_t *array, fc_json_value_t *value)
{
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON value.");
    }
    fc_json_check(json_array_append_value(array, value));
}

void fc_json_array_append_string(fc_json_array_t *array, const char *value)
{
    fc_json_check(json_array_append_string(array, value ? value : ""));
}

void fc_json_array_append_number(fc_json_array_t *array, double value)
{
    fc_json_check(json_array_append_number(array, fc_json_number(value)));
}

fc_json_array_t *fc_json_array_append_array(fc_json_array_t *array)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_array_append_value(array, value);
    return fc_json_get_array(value);
}

fc_json_value_t *fc_json_number_array_value(const double *values, size_t n)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_array_t *array = fc_json_get_array(value);
    for (size_t i = 0; i < n; i++) {
        fc_json_array_append_number(array, values ? values[i] : 0.0);
    }
    return value;
}

void fc_json_entries_append_value(fc_json_array_t *entries, const char *key, fc_json_value_t *value)
{
    fc_json_array_t *entry = fc_json_array_append_array(entries);
    fc_json_array_append_string(entry, key);
    fc_json_array_append_value(entry, value);
}
