#include "core/fc_time.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static bool parse_decimal(const char *text, size_t offset, size_t digits, int *out)
{
    int value = 0;
    for (size_t i = 0; i < digits; ++i) {
        char character = text[offset + i];
        if (character < '0' || character > '9') {
            return false;
        }
        value = value * 10 + (character - '0');
    }
    *out = value;
    return true;
}

static bool is_leap_year(int year)
{
    return (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);
}

static bool valid_date(int year, int month, int day)
{
    static const unsigned char days_per_month[] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
    if (month < 1 || month > 12 || day < 1) {
        return false;
    }
    int days = days_per_month[month - 1];
    if (month == 2 && is_leap_year(year)) {
        ++days;
    }
    return day <= days;
}

// Days since 1970-01-01 in the proleptic Gregorian calendar.  This is the
// March-based 400-year-cycle calculation used by timegm implementations, but
// only for already validated calendar fields.
static int64_t days_from_civil(int year, int month, int day)
{
    int64_t adjusted_year = year;
    if (month <= 2) {
        --adjusted_year;
    }
    int64_t era = adjusted_year >= 0 ? adjusted_year / 400 : (adjusted_year - 399) / 400;
    unsigned year_of_era = (unsigned)(adjusted_year - era * 400);
    unsigned march_month = (unsigned)(month > 2 ? month - 3 : month + 9);
    unsigned day_of_year = (153 * march_month + 2) / 5 + (unsigned)day - 1;
    unsigned day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    return era * 146097 + (int64_t)day_of_era - 719468;
}

static bool parse_iso_timestamp(const char *text, double *out)
{
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    int second = 0;
    if (!parse_decimal(text, 0, 4, &year)
        || !parse_decimal(text, 5, 2, &month)
        || !parse_decimal(text, 8, 2, &day)
        || !parse_decimal(text, 11, 2, &hour)
        || !parse_decimal(text, 14, 2, &minute)
        || !parse_decimal(text, 17, 2, &second)
        || !valid_date(year, month, day)
        || hour > 23 || minute > 59 || second > 59) {
        return false;
    }

    int64_t days = days_from_civil(year, month, day);
    int64_t seconds = days * 86400 + hour * 3600 + minute * 60 + second;
    double timestamp = (double)seconds;
    if (!isfinite(timestamp)) {
        return false;
    }
    *out = timestamp;
    return true;
}

bool fc_time_parse_timestamp(const char *text, double *out)
{
    if (NULL == text || NULL == out || text[0] == '\0') {
        return false;
    }

    size_t len = strlen(text);
    if (len >= 19 && text[4] == '-' && text[7] == '-' && text[10] == 'T'
        && text[13] == ':' && text[16] == ':') {
        return parse_iso_timestamp(text, out);
    }

    errno = 0;
    char *end = NULL;
    double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value)) {
        return false;
    }
    *out = value;
    return true;
}
