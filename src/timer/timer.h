#ifndef TIMER_H
#define TIMER_H

#include <windows.h>

typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    LARGE_INTEGER frequency;
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
