#ifndef FASTCAT_ARGP_H
#define FASTCAT_ARGP_H

/*
 * MinGW's default __printf__ format checker follows the legacy Microsoft
 * printf convention, under which the C99 %z length modifier is unsupported.
 * argp-standalone uses that checker for argp_error and related functions.
 * UCRT itself supports %zu, so select GCC's C99-compatible checker before
 * including the vendored argp header.
 */
#if defined(__MINGW32__)
#define PRINTF_STYLE(f, a) __attribute__((__format__(gnu_printf, f, a)))
#endif

#include <argp.h>

#endif
