#ifndef TIMER_H
#define TIMER_H

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

typedef struct {
#ifdef _WIN32
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    LARGE_INTEGER frequency;
#else
    struct timespec start;
    struct timespec end;
#endif
} Timer;

void timer_start(Timer* t);

void timer_stop(Timer* t);

/* Вернуть прошедшее время в секундах. */
double timer_elapsed_s(const Timer* t);

/* Вернуть прошедшее время в миллисекундах. */
double timer_elapsed_ms(const Timer* t);

/* Вернуть прошедшее время в микросекундах. */
double timer_elapsed_us(const Timer* t);

#endif
