#include "fc_log.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    fc_log_level level;
    char module[11]; // 10 chars + NUL
    FILE *target;
    fc_log_color_mode color;
    pthread_mutex_t lock;
} fc_log_state;

static fc_log_state g_log = {
    .level = FC_LOG_INFO,
    .module = "fastcat",
    .target = NULL,
    .color = FC_LOG_COLOR_AUTO,
    .lock = PTHREAD_MUTEX_INITIALIZER
};

static const char *LEVEL_NAME[] = { "ERROR", "WARN ", "INFO ", "DEBUG" };
#define LEVEL_COUNT 4

static void build_module(const char *src, char out[11])
{
    memset(out, ' ', 10);
    out[10] = '\0';
    if (!src) {
        memcpy(out, g_log.module, 10);
        return;
    }
    size_t len = strlen(src);
    if (len > 10) {
        len = 10;
    }
    memcpy(out, src, len);
}

fc_log_level fc_log_level_from_string(const char *s, fc_log_level default_level)
{
    if (!s) {
        return default_level;
    }
    if (strcasecmp(s, "error") == 0) {
        return FC_LOG_ERROR;
    }
    if (strcasecmp(s, "warn") == 0 || strcasecmp(s, "warning") == 0) {
        return FC_LOG_WARN;
    }
    if (strcasecmp(s, "info") == 0) {
        return FC_LOG_INFO;
    }
    if (strcasecmp(s, "debug") == 0) {
        return FC_LOG_DEBUG;
    }
    return default_level;
}

void fc_log_init(const fc_log_config *cfg)
{
    fc_log_config def = {
        .level = FC_LOG_INFO,
        .module = "fastcat",
        .target = NULL,
        .color = FC_LOG_COLOR_AUTO
    };
    if (!cfg) {
        cfg = &def;
    }

    g_log.level = cfg->level;

    build_module(cfg->module, g_log.module);

    g_log.target = cfg->target ? fopen(cfg->target, "w") : stderr;
    if (!g_log.target) {
        g_log.target = stderr;
    }

    g_log.color = cfg->color;
    if (g_log.color == FC_LOG_COLOR_AUTO) {
        int tty = g_log.target ? isatty(fileno(g_log.target)) : 0;
        g_log.color = cfg->target ? FC_LOG_COLOR_OFF : (tty ? FC_LOG_COLOR_ON : FC_LOG_COLOR_OFF);
    }
}

void fc_log_set_level(fc_log_level level)
{
    g_log.level = level;
}

void fc_log_destroy(void)
{
    pthread_mutex_lock(&g_log.lock);
    if (g_log.target && g_log.target != stderr) {
        fclose(g_log.target);
    }
    g_log.target = NULL;
    pthread_mutex_unlock(&g_log.lock);
}

static void format_time(char *buf, size_t n)
{
    struct timespec ts;
    struct tm tm;
    clock_gettime(CLOCK_REALTIME, &ts);
    localtime_r(&ts.tv_sec, &tm);
    snprintf(buf, n, "%02d:%02d:%02d", tm.tm_hour, tm.tm_min, tm.tm_sec);
}

void fc_logv(fc_log_level level, const char *fmt, va_list ap)
{
    fc_logmv(level, NULL, fmt, ap);
}

void fc_logmv(fc_log_level level, const char *module, const char *fmt, va_list ap)
{
    if (level > g_log.level) {
        return;
    }
    if (g_log.target == NULL) {
        g_log.target = stderr;
    }
    pthread_mutex_lock(&g_log.lock);

    char tbuf[16];
    format_time(tbuf, sizeof(tbuf));

    char modbuf[11];
    build_module(module, modbuf);

    const char *lvl = (level >= 0 && level < (int)LEVEL_COUNT) ? LEVEL_NAME[level] : "?????";
    const char *c_start = "";
    const char *c_end = "";
    if (g_log.color == FC_LOG_COLOR_ON) {
        switch (level) {
        case FC_LOG_ERROR:
            c_start = "\x1b[1;31m";
            break;
        case FC_LOG_WARN:
            c_start = "\x1b[1;33m";
            break;
        case FC_LOG_INFO:
            c_start = "\x1b[1;32m";
            break;
        case FC_LOG_DEBUG:
            c_start = "\x1b[1;35m";
            break;
        }
        c_end = "\x1b[0m";
    }

    fprintf(g_log.target, "[%s - %s] %s%s%s ", tbuf, modbuf, c_start, lvl, c_end);
    vfprintf(g_log.target, fmt, ap);
    fputc('\n', g_log.target);
    fflush(g_log.target);

    pthread_mutex_unlock(&g_log.lock);
}

void fc_log(fc_log_level level, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fc_logv(level, fmt, ap);
    va_end(ap);
}

void fc_logm(fc_log_level level, const char *module, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fc_logmv(level, module, fmt, ap);
    va_end(ap);
}
