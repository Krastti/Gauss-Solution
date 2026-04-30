#include "timer.h"
#include "../logging/logger.h"

#ifdef _WIN32
static double timer_ticks_elapsed_s(const Timer* t) {
    double ticks;

    if (t->frequency.QuadPart == 0) {
        LOG_ERROR("%s", "Частота таймера равна нулю!");
        return -1.0;
    }

    ticks = (double)(t->end.QuadPart - t->start.QuadPart);
    return ticks / (double)t->frequency.QuadPart;
}
#else
static double timer_timespec_elapsed_s(const Timer* t) {
    const time_t seconds = t->end.tv_sec - t->start.tv_sec;
    const long nanoseconds = t->end.tv_nsec - t->start.tv_nsec;

    return (double)seconds + (double)nanoseconds / 1000000000.0;
}
#endif

void timer_start(Timer* t) {
    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return;
    }

#ifdef _WIN32
    QueryPerformanceFrequency(&t->frequency);
    QueryPerformanceCounter(&t->start);
#else
    clock_gettime(CLOCK_MONOTONIC, &t->start);
#endif
}

void timer_stop(Timer* t) {
    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return;
    }

#ifdef _WIN32
    QueryPerformanceCounter(&t->end);
#else
    clock_gettime(CLOCK_MONOTONIC, &t->end);
#endif
}

double timer_elapsed_s(const Timer* t) {
    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return -1.0;
    }

#ifdef _WIN32
    return timer_ticks_elapsed_s(t);
#else
    return timer_timespec_elapsed_s(t);
#endif
}

double timer_elapsed_ms(const Timer* t) {
    return timer_elapsed_s(t) * 1000.0;
}

double timer_elapsed_us(const Timer* t) {
    return timer_elapsed_s(t) * 1000000.0;
}
