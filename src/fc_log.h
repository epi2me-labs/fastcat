#ifndef FC_LOG_H
#define FC_LOG_H

#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>

typedef enum {
    FC_LOG_ERROR = 0,
    FC_LOG_WARN  = 1,
    FC_LOG_INFO  = 2,
    FC_LOG_DEBUG = 3
} fc_log_level;

typedef enum {
    FC_LOG_COLOR_OFF = 0,
    FC_LOG_COLOR_ON  = 1,
    FC_LOG_COLOR_AUTO = 2
} fc_log_color_mode;

typedef struct {
    fc_log_level level;          // default FC_LOG_INFO
    const char *module;          // padded/truncated to 10 chars
    const char *target;          // NULL => stderr
    fc_log_color_mode color;     // on/off/auto
} fc_log_config;

// Initialize logger once early in main.
void fc_log_init(const fc_log_config *cfg);
void fc_log_set_level(fc_log_level level);
void fc_log_destroy(void);
fc_log_level fc_log_level_from_string(const char *s, fc_log_level default_level);
void fc_logv(fc_log_level level, const char *fmt, va_list ap);
void fc_log(fc_log_level level, const char *fmt, ...);
void fc_logm(fc_log_level level, const char *module, const char *fmt, ...);
void fc_logmv(fc_log_level level, const char *module, const char *fmt, va_list ap);

#define FC_LOG_ERR(...)   fc_log(FC_LOG_ERROR, __VA_ARGS__)
#define FC_LOG_WARN(...)  fc_log(FC_LOG_WARN,  __VA_ARGS__)
#define FC_LOG_INFO(...)  fc_log(FC_LOG_INFO,  __VA_ARGS__)
#define FC_LOG_DEBUG(...) fc_log(FC_LOG_DEBUG, __VA_ARGS__)
#define FC_LOG_ERROR(...) FC_LOG_ERR(__VA_ARGS__)

#define FC_LOGM_ERR(mod, ...)   fc_logm(FC_LOG_ERROR, mod, __VA_ARGS__)
#define FC_LOGM_WARN(mod, ...)  fc_logm(FC_LOG_WARN,  mod, __VA_ARGS__)
#define FC_LOGM_INFO(mod, ...)  fc_logm(FC_LOG_INFO,  mod, __VA_ARGS__)
#define FC_LOGM_DEBUG(mod, ...) fc_logm(FC_LOG_DEBUG, mod, __VA_ARGS__)
#define FC_LOGM_ERROR(mod, ...) FC_LOGM_ERR(mod, __VA_ARGS__)

#endif /* FC_LOG_H */
