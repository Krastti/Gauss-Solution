#include "logger.h"

#include <stdarg.h>
#include <time.h>

static LogLevel current_level = LOG_DEBUG;

static const char* LEVEL_NAMES[] = {
    "DEBUG",
    "INFO",
    "WARN",
    "ERROR"
};

LogLevel logger_current_level(void) {
  return current_level;
}

void logger_set_level(LogLevel level) {
    if (level < LOG_DEBUG || level > LOG_ERROR) {
        printf("[LOGGER] Некорректный уровень логирования: %d\n", level);
        return;
    }

    current_level = level;
    printf("[LOGGER] Уровень логирования установлен: %s\n", LEVEL_NAMES[level]);
}

void logger_log(LogLevel level, const char* fmt, ...) {
    time_t now;
    struct tm* t;
    char timebuf[20];
    va_list args;

    if (level < LOG_DEBUG || level > LOG_ERROR) {
        printf("[LOGGER] Некорректный уровень сообщения: %d\n", level);
        return;
    }

    if (level < current_level) {
        return;
    }

    now = time(NULL);
    t = localtime(&now);
    strftime(timebuf, sizeof(timebuf), "%H:%M:%S", t);

    printf("[%s] [%s] ", timebuf, LEVEL_NAMES[level]);

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    printf("\n");
}
