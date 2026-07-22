#include <math.h>
#include <stdint.h>
#include <string.h>

#include "unity.h"
#include "main.h"
#include "core/fc_json.h"
#include "vendor/parson/parson.h"


static JSON_Value *round_trip_json_value(fc_json_value_t *value)
{
    char *text = fc_json_serialize_to_string(value);
    TEST_ASSERT_NOT_NULL(text);
    JSON_Value *parsed = json_parse_string(text);
    json_free_serialized_string(text);
    fc_json_value_free(value);
    TEST_ASSERT_NOT_NULL(parsed);
    return parsed;
}


static void assert_shape(const JSON_Object *object, size_t rows, size_t cols)
{
    JSON_Array *shape = json_object_get_array(object, "shape");
    TEST_ASSERT_NOT_NULL(shape);
    TEST_ASSERT_EQUAL_size_t(2, json_array_get_count(shape));
    TEST_ASSERT_EQUAL_UINT64(rows, (uint64_t)json_array_get_number(shape, 0));
    TEST_ASSERT_EQUAL_UINT64(cols, (uint64_t)json_array_get_number(shape, 1));
}


static void assert_adaptive_column_dtype(size_t max_value, const char *expected_dtype, size_t width)
{
    JSON_Value *parsed = round_trip_json_value(fc_json_column_uint_adaptive("value", &max_value, 1));
    JSON_Object *object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);

    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    TEST_ASSERT_EQUAL_STRING("value", json_object_get_string(object, "name"));
    TEST_ASSERT_EQUAL_STRING(expected_dtype, json_object_get_string(object, "dtype"));
    TEST_ASSERT_FALSE(json_object_has_value(object, "shape"));
    TEST_ASSERT_EQUAL_size_t(4 * ((width + 2) / 3), strlen(json_object_get_string(object, "__ndarray__")));

    json_value_free(parsed);
}


FC_TEST(json_suite, test_adaptive_uint_column_uses_smallest_dtype)
{
    assert_adaptive_column_dtype(UINT8_MAX, "uint8", sizeof(uint8_t));
    assert_adaptive_column_dtype((size_t)UINT8_MAX + 1, "uint16", sizeof(uint16_t));
    assert_adaptive_column_dtype((size_t)UINT16_MAX + 1, "uint32", sizeof(uint32_t));
    if (SIZE_MAX > UINT32_MAX) {
        assert_adaptive_column_dtype((size_t)UINT32_MAX + 1, "uint64", sizeof(uint64_t));
    }
}


FC_TEST(json_suite, test_dataframe_columns_omit_shape_and_preserve_name)
{
    const size_t counts[] = { 1, 255 };
    JSON_Value *parsed = round_trip_json_value(fc_json_column_uint_adaptive("read_count", counts, 2));
    JSON_Object *object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);

    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    TEST_ASSERT_EQUAL_STRING("read_count", json_object_get_string(object, "name"));
    TEST_ASSERT_EQUAL_STRING("uint8", json_object_get_string(object, "dtype"));
    TEST_ASSERT_FALSE(json_object_has_value(object, "shape"));
    TEST_ASSERT_EQUAL_STRING("Af8=", json_object_get_string(object, "__ndarray__"));
    json_value_free(parsed);

    const double means[] = { 1.5, 2.5 };
    parsed = round_trip_json_value(fc_json_column_f64("length_mean", means, 2));
    object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    TEST_ASSERT_EQUAL_STRING("length_mean", json_object_get_string(object, "name"));
    TEST_ASSERT_EQUAL_STRING("float64", json_object_get_string(object, "dtype"));
    TEST_ASSERT_FALSE(json_object_has_value(object, "shape"));
    TEST_ASSERT_EQUAL_size_t(24, strlen(json_object_get_string(object, "__ndarray__")));
    json_value_free(parsed);
}


FC_TEST(json_suite, test_matrix_ndarrays_include_shape_dtype_and_correct_payload_size)
{
    const double floats[] = { 1.5, 2.5, 3.5, 4.5 };
    JSON_Value *parsed = round_trip_json_value(fc_json_ndarray_f64(floats, 2, 2));
    JSON_Object *object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    assert_shape(object, 2, 2);
    TEST_ASSERT_EQUAL_STRING("float64", json_object_get_string(object, "dtype"));
    TEST_ASSERT_EQUAL_size_t(44, strlen(json_object_get_string(object, "__ndarray__")));
    json_value_free(parsed);

    const size_t raw_counts[] = { 1, 2, 3, 4 };
    parsed = round_trip_json_value(fc_json_ndarray_uint_adaptive(raw_counts, 2, 2));
    object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    assert_shape(object, 2, 2);
    TEST_ASSERT_EQUAL_STRING("uint8", json_object_get_string(object, "dtype"));
    TEST_ASSERT_EQUAL_size_t(8, strlen(json_object_get_string(object, "__ndarray__")));
    json_value_free(parsed);

    const double counts[] = { -1.0, NAN, 1.2, 300.0 };
    parsed = round_trip_json_value(fc_json_ndarray_f64_to_uint_adaptive(counts, 2, 2));
    object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    assert_shape(object, 2, 2);
    TEST_ASSERT_EQUAL_STRING("uint16", json_object_get_string(object, "dtype"));
    TEST_ASSERT_EQUAL_size_t(12, strlen(json_object_get_string(object, "__ndarray__")));
    json_value_free(parsed);

    parsed = round_trip_json_value(fc_json_ndarray_f64_to_uint_adaptive(NULL, 0, 3));
    object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    assert_shape(object, 0, 3);
    TEST_ASSERT_EQUAL_STRING("uint8", json_object_get_string(object, "dtype"));
    TEST_ASSERT_EQUAL_STRING("", json_object_get_string(object, "__ndarray__"));
    json_value_free(parsed);
}


FC_TEST(json_suite, test_empty_adaptive_column_has_uint8_empty_payload)
{
    JSON_Value *parsed = round_trip_json_value(fc_json_column_uint_adaptive("count", NULL, 0));
    JSON_Object *object = json_value_get_object(parsed);
    TEST_ASSERT_NOT_NULL(object);
    TEST_ASSERT_EQUAL_size_t(3, json_object_get_count(object));
    TEST_ASSERT_EQUAL_STRING("count", json_object_get_string(object, "name"));
    TEST_ASSERT_EQUAL_STRING("uint8", json_object_get_string(object, "dtype"));
    TEST_ASSERT_FALSE(json_object_has_value(object, "shape"));
    TEST_ASSERT_EQUAL_STRING("", json_object_get_string(object, "__ndarray__"));
    json_value_free(parsed);
}
