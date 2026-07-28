#include "core/fc_json.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "htslib/kstring.h"

#include "core/fc_common.h"
#include "core/fc_log.h"
#include "vendor/parson/parson.h"

/* Allocate Parson-owned memory, treating a zero-byte request as one byte. */
static void *fc_json_malloc(size_t size)
{
    void *ptr = malloc(size ? size : 1);
    if (!ptr) {
        FC_LOG_FATAL("Failed to allocate JSON memory.");
    }
    return ptr;
}

/* Release memory allocated through Parson's configured allocator. */
static void fc_json_free(void *ptr)
{
    free(ptr);
}

/* Configure Parson once so its allocation failures use fastcat's fatal path. */
static void fc_json_init(void)
{
    static bool initialized = false;
    if (!initialized) {
        json_set_allocation_functions(fc_json_malloc, fc_json_free);
        initialized = true;
    }
}

/* Turn a failed Parson mutation into a fatal fastcat error. */
static void fc_json_check(JSON_Status status)
{
    if (status != JSONSuccess) {
        FC_LOG_FATAL("Failed to build JSON.");
    }
}

/* JSON has no NaN or infinity literals, so replace non-finite values with zero. */
static double fc_json_number(double value)
{
    return isfinite(value) ? value : 0.0;
}

/* Append RFC 4648 base64 for an arbitrary binary payload. */
static void fc_json_append_base64(kstring_t *ks, const uint8_t *data, size_t n)
{
    static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    for (size_t i = 0; i < n; i += 3) {
        uint32_t a = data[i];
        uint32_t b = (i + 1 < n) ? data[i + 1] : 0;
        uint32_t c = (i + 2 < n) ? data[i + 2] : 0;
        uint32_t triple = (a << 16) | (b << 8) | c;

        kputc(table[(triple >> 18) & 0x3f], ks);
        kputc(table[(triple >> 12) & 0x3f], ks);
        kputc((i + 1 < n) ? table[(triple >> 6) & 0x3f] : '=', ks);
        kputc((i + 2 < n) ? table[triple & 0x3f] : '=', ks);
    }
}

/* Encode bytes as an owned, NUL-terminated base64 string. */
static char *fc_json_bytes_base64(const uint8_t *bytes, size_t n)
{
    kstring_t b64 = { 0, 0, NULL };
    fc_json_append_base64(&b64, bytes, n);
    return ks_release(&b64);
}

/* Create an empty JSON object value. */
fc_json_value_t *fc_json_object_value(void)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_object();
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON object.");
    }
    return value;
}

/* Create an empty JSON array value. */
fc_json_value_t *fc_json_array_value(void)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_array();
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON array.");
    }
    return value;
}

/* Create a JSON string, mapping NULL input to the empty string. */
fc_json_value_t *fc_json_string_value(const char *text)
{
    fc_json_init();
    fc_json_value_t *value = json_value_init_string(text ? text : "");
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON string.");
    }
    return value;
}

/* Create a finite JSON number value. */
fc_json_value_t *fc_json_number_value(double value)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_number(fc_json_number(value));
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON number.");
    }
    return json_value;
}

/* Create a JSON boolean value. */
fc_json_value_t *fc_json_bool_value(int value)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_boolean(value);
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON boolean.");
    }
    return json_value;
}

/* Create a JSON null value. */
fc_json_value_t *fc_json_null_value(void)
{
    fc_json_init();
    fc_json_value_t *json_value = json_value_init_null();
    if (!json_value) {
        FC_LOG_FATAL("Failed to allocate JSON null.");
    }
    return json_value;
}

/* Parse one JSON document, failing fast on invalid input. */
fc_json_value_t *fc_json_parse_string(const char *json)
{
    fc_json_init();
    fc_json_value_t *value = json_parse_string(json ? json : "");
    if (!value) {
        FC_LOG_FATAL("Failed to parse JSON.");
    }
    return value;
}

/* Return an object's Parson handle, failing if value is not an object. */
fc_json_object_t *fc_json_get_object(fc_json_value_t *value)
{
    fc_json_object_t *object = json_object(value);
    if (!object) {
        FC_LOG_FATAL("Failed to access JSON object.");
    }
    return object;
}

/* Return an array's Parson handle, failing if value is not an array. */
fc_json_array_t *fc_json_get_array(fc_json_value_t *value)
{
    fc_json_array_t *array = json_array(value);
    if (!array) {
        FC_LOG_FATAL("Failed to access JSON array.");
    }
    return array;
}

/* Free a JSON value and every value it owns. */
void fc_json_value_free(fc_json_value_t *value)
{
    json_value_free(value);
}

/* Serialize a JSON tree to an owned compact string. */
char *fc_json_serialize_to_string(const fc_json_value_t *value)
{
    fc_json_init();
    char *serialized = json_serialize_to_string(value);
    if (!serialized) {
        FC_LOG_FATAL("Failed to serialise JSON.");
    }
    return serialized;
}

/* Transfer ownership of value into object under key; NULL writes JSON null. */
void fc_json_object_set_value(fc_json_object_t *object, const char *key, fc_json_value_t *value)
{
    if (!value) {
        fc_json_object_set_null(object, key);
        return;
    }
    fc_json_check(json_object_set_value(object, key, value));
}

/* Store a string member, mapping NULL input to the empty string. */
void fc_json_object_set_string(fc_json_object_t *object, const char *key, const char *value)
{
    fc_json_check(json_object_set_string(object, key, value ? value : ""));
}

/* Store a finite numeric member. */
void fc_json_object_set_number(fc_json_object_t *object, const char *key, double value)
{
    fc_json_check(json_object_set_number(object, key, fc_json_number(value)));
}

/* Store a null member. */
void fc_json_object_set_null(fc_json_object_t *object, const char *key)
{
    fc_json_check(json_object_set_null(object, key));
}

/* Create, attach, and return a child object. */
fc_json_object_t *fc_json_object_set_object(fc_json_object_t *object, const char *key)
{
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_set_value(object, key, value);
    return fc_json_get_object(value);
}

/* Create, attach, and return a child array. */
fc_json_array_t *fc_json_object_set_array(fc_json_object_t *object, const char *key)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_object_set_value(object, key, value);
    return fc_json_get_array(value);
}

/* Transfer ownership of value to the end of array. */
void fc_json_array_append_value(fc_json_array_t *array, fc_json_value_t *value)
{
    if (!value) {
        FC_LOG_FATAL("Failed to allocate JSON value.");
    }
    fc_json_check(json_array_append_value(array, value));
}

/* Append a string, mapping NULL input to the empty string. */
void fc_json_array_append_string(fc_json_array_t *array, const char *value)
{
    fc_json_check(json_array_append_string(array, value ? value : ""));
}

/* Append a finite JSON number. */
void fc_json_array_append_number(fc_json_array_t *array, double value)
{
    fc_json_check(json_array_append_number(array, fc_json_number(value)));
}

/* Create, append, and return a nested array. */
fc_json_array_t *fc_json_array_append_array(fc_json_array_t *array)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_array_append_value(array, value);
    return fc_json_get_array(value);
}

/* Encode n doubles as a regular JSON number array. */
fc_json_value_t *fc_json_number_array_value(const double *values, size_t n)
{
    fc_json_value_t *value = fc_json_array_value();
    fc_json_array_t *array = fc_json_get_array(value);
    for (size_t i = 0; i < n; i++) {
        fc_json_array_append_number(array, values ? values[i] : 0.0);
    }
    return value;
}

/* Build a self-describing little-endian binary array. */
static fc_json_value_t *fc_json_ndarray_value(
    const char *name, const char *dtype, const uint8_t *bytes, size_t n_bytes, size_t element_size,
    const size_t *shape, size_t n_dims)
{
    const uint8_t *encoded_bytes = bytes;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    uint8_t *little_endian_bytes = n_bytes ? fc_alloc(n_bytes, sizeof(*little_endian_bytes), "JSON little-endian array") : NULL;
    for (size_t offset = 0; offset < n_bytes; offset += element_size) {
        for (size_t byte = 0; byte < element_size; ++byte) {
            little_endian_bytes[offset + byte] = bytes[offset + element_size - byte - 1];
        }
    }
    encoded_bytes = little_endian_bytes;
#else
    (void)element_size;
#endif
    char *b64 = fc_json_bytes_base64(encoded_bytes, n_bytes);
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);

    if (name) {
        fc_json_object_set_string(object, "name", name);
    }
    fc_json_object_set_string(object, "kind", "array");
    if (n_dims > 0) {
        fc_json_array_t *json_shape = fc_json_object_set_array(object, "shape");
        for (size_t i = 0; i < n_dims; i++) {
            fc_json_array_append_number(json_shape, (double)shape[i]);
        }
    }
    fc_json_object_set_string(object, "dtype", dtype);
    fc_json_object_set_string(object, "data", b64);

    free(b64);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    free(little_endian_bytes);
#endif
    return value;
}

/* Mark an object as a column-oriented dataframe and return its columns array. */
fc_json_array_t *fc_json_object_set_dataframe(
    fc_json_object_t *object, size_t rows, size_t columns)
{
    fc_json_object_set_string(object, "kind", "dataframe");
    fc_json_array_t *shape = fc_json_object_set_array(object, "shape");
    fc_json_array_append_number(shape, (double)rows);
    fc_json_array_append_number(shape, (double)columns);
    return fc_json_object_set_array(object, "columns");
}

/* Return the element count for a two-dimensional ndarray, rejecting overflow. */
static size_t fc_json_ndarray_2d_count(size_t rows, size_t cols)
{
    if (rows != 0 && cols > SIZE_MAX / rows) {
        FC_LOG_FATAL("JSON ndarray is too large.");
    }
    return rows * cols;
}

/* Pack size_t values into the narrowest unsigned dtype, then build an ndarray object. */
static fc_json_value_t *fc_json_adaptive_uint_value(
    const char *name, const size_t *values, size_t n, const size_t *shape, size_t n_dims)
{
    if (n == 0) {
        return fc_json_ndarray_value(name, "<u1", NULL, 0, sizeof(uint8_t), shape, n_dims);
    }

    size_t max_value = 0;
    for (size_t i = 0; i < n; i++) {
        if (values[i] > max_value) {
            max_value = values[i];
        }
    }

    const char *dtype;
    size_t element_size;
    if (max_value <= UINT8_MAX) {
        dtype = "<u1";
        element_size = sizeof(uint8_t);
    } else if (max_value <= UINT16_MAX) {
        dtype = "<u2";
        element_size = sizeof(uint16_t);
    } else if (max_value <= UINT32_MAX) {
        dtype = "<u4";
        element_size = sizeof(uint32_t);
    } else {
        dtype = "<u8";
        element_size = sizeof(uint64_t);
    }

    void *packed = fc_alloc(n, element_size, "adaptive JSON ndarray");
    if (element_size == sizeof(uint8_t)) {
        for (size_t i = 0; i < n; i++) {
            ((uint8_t *)packed)[i] = (uint8_t)values[i];
        }
    } else if (element_size == sizeof(uint16_t)) {
        for (size_t i = 0; i < n; i++) {
            ((uint16_t *)packed)[i] = (uint16_t)values[i];
        }
    } else if (element_size == sizeof(uint32_t)) {
        for (size_t i = 0; i < n; i++) {
            ((uint32_t *)packed)[i] = (uint32_t)values[i];
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            ((uint64_t *)packed)[i] = (uint64_t)values[i];
        }
    }

    fc_json_value_t *value = fc_json_ndarray_value(
        name, dtype, (const uint8_t *)packed, n * element_size, element_size, shape, n_dims);
    free(packed);
    return value;
}

/* Encode a standalone two-dimensional ndarray using the narrowest unsigned dtype. */
fc_json_value_t *fc_json_ndarray_uint_adaptive(const size_t *values, size_t rows, size_t cols)
{
    size_t n = fc_json_ndarray_2d_count(rows, cols);
    size_t shape[] = { rows, cols };
    return fc_json_adaptive_uint_value(NULL, values, n, shape, 2);
}

/* Encode a standalone two-dimensional little-endian float64 ndarray. */
fc_json_value_t *fc_json_ndarray_f64(const double *values, size_t rows, size_t cols)
{
    size_t n = fc_json_ndarray_2d_count(rows, cols);
    if (n > SIZE_MAX / sizeof(double)) {
        FC_LOG_FATAL("JSON ndarray is too large.");
    }

    size_t shape[] = { rows, cols };
    return fc_json_ndarray_value(
        NULL, "<f8", (const uint8_t *)values, n * sizeof(*values), sizeof(*values), shape, 2);
}

/* Encode a two-dimensional float64 ndarray using Bokeh 3's JSON wire format. */
fc_json_value_t *fc_json_bokeh_ndarray_f64(const double *values, size_t rows, size_t cols)
{
    size_t n = fc_json_ndarray_2d_count(rows, cols);
    if (n > SIZE_MAX / sizeof(double)) {
        FC_LOG_FATAL("Bokeh JSON ndarray is too large.");
    }

    char *b64 = fc_json_bytes_base64((const uint8_t *)values, n * sizeof(*values));
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);
    fc_json_object_set_string(object, "type", "ndarray");
    fc_json_object_t *array = fc_json_object_set_object(object, "array");
    fc_json_object_set_string(array, "type", "bytes");
    fc_json_object_set_string(array, "data", b64);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    fc_json_object_set_string(object, "order", "big");
#else
    fc_json_object_set_string(object, "order", "little");
#endif
    fc_json_object_set_string(object, "dtype", "float64");
    fc_json_array_t *shape = fc_json_object_set_array(object, "shape");
    fc_json_array_append_number(shape, (double)rows);
    fc_json_array_append_number(shape, (double)cols);
    free(b64);
    return value;
}

/* Encode a two-dimensional float32 ndarray using Bokeh 3's JSON wire format. */
fc_json_value_t *fc_json_bokeh_ndarray_f32(const double *values, size_t rows, size_t cols)
{
    size_t n = fc_json_ndarray_2d_count(rows, cols);
    if (n > SIZE_MAX / sizeof(float)) {
        FC_LOG_FATAL("Bokeh JSON ndarray is too large.");
    }

    float *packed = fc_alloc(n, sizeof(*packed), "Bokeh float32 ndarray");
    for (size_t i = 0; i < n; ++i) {
        packed[i] = (float)values[i];
    }

    char *b64 = fc_json_bytes_base64((const uint8_t *)packed, n * sizeof(*packed));
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);
    fc_json_object_set_string(object, "type", "ndarray");
    fc_json_object_t *array = fc_json_object_set_object(object, "array");
    fc_json_object_set_string(array, "type", "bytes");
    fc_json_object_set_string(array, "data", b64);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    fc_json_object_set_string(object, "order", "big");
#else
    fc_json_object_set_string(object, "order", "little");
#endif
    fc_json_object_set_string(object, "dtype", "float32");
    fc_json_array_t *shape = fc_json_object_set_array(object, "shape");
    fc_json_array_append_number(shape, (double)rows);
    fc_json_array_append_number(shape, (double)cols);
    free(b64);
    free(packed);
    return value;
}

/* Encode a named adaptive unsigned dataframe column. */
fc_json_value_t *fc_json_column_uint_adaptive(const char *name, const size_t *values, size_t n)
{
    size_t shape[] = { n };
    return fc_json_adaptive_uint_value(name, values, n, shape, 1);
}

/* Encode a named float64 dataframe column. */
fc_json_value_t *fc_json_column_f64(const char *name, const double *values, size_t n)
{
    size_t shape[] = { n };
    return fc_json_ndarray_value(
        name, "<f8", (const uint8_t *)values, n * sizeof(*values), sizeof(*values), shape, 1);
}

/* Encode a named dataframe string column directly, without binary ndarray packing. */
fc_json_value_t *fc_json_column_string(const char *name, const char *const *values, size_t n)
{
    fc_json_value_t *value = fc_json_object_value();
    fc_json_object_t *object = fc_json_get_object(value);
    fc_json_object_set_string(object, "name", name);
    fc_json_object_set_string(object, "kind", "array");
    fc_json_object_set_string(object, "dtype", "string");
    fc_json_array_t *shape = fc_json_object_set_array(object, "shape");
    fc_json_array_append_number(shape, (double)n);
    fc_json_array_t *data = fc_json_object_set_array(object, "data");
    for (size_t i = 0; i < n; ++i) {
        fc_json_array_append_string(data, values ? values[i] : "");
    }
    return value;
}

/* Quantize a double matrix into the smallest unsigned dtype and attach its two-dimensional shape. */
fc_json_value_t *fc_json_ndarray_f64_to_uint_adaptive(const double *values, size_t rows, size_t cols)
{
    size_t n = fc_json_ndarray_2d_count(rows, cols);
    if (n == 0) {
        return fc_json_ndarray_uint_adaptive(NULL, rows, cols);
    }

    size_t *counts = fc_alloc(n, sizeof(*counts), "quantized heatmap counts");
    for (size_t i = 0; i < n; i++) {
        double value = values[i];
        if (!isfinite(value) || value < 0.0) {
            counts[i] = 0;
        } else if (value > UINT32_MAX) {
            counts[i] = UINT32_MAX;
        } else {
            counts[i] = (size_t)llround(value);
        }
    }

    fc_json_value_t *result = fc_json_ndarray_uint_adaptive(counts, rows, cols);
    free(counts);
    return result;
}

/* Append a [key, value] pair to a JSON entries array. */
void fc_json_entries_append_value(fc_json_array_t *entries, const char *key, fc_json_value_t *value)
{
    fc_json_array_t *entry = fc_json_array_append_array(entries);
    fc_json_array_append_string(entry, key);
    fc_json_array_append_value(entry, value);
}
