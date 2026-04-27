#include "timer.h"
#include "../logging/logger.h"

void timer_start(Timer* t) {
    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return;
    }

    QueryPerformanceFrequency(&t->frequency);
    QueryPerformanceCounter(&t->start);
}

void timer_stop(Timer* t) {
    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return;
    }

    QueryPerformanceCounter(&t->end);
}

double timer_elapsed_s(const Timer* t) {
    double ticks;

    if (t == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель на таймер!");
        return -1.0;
    }

    if (t->frequency.QuadPart == 0) {
        LOG_ERROR("%s", "Частота таймера равна нулю!");
        return -1.0;
    }

    ticks = (double)(t->end.QuadPart - t->start.QuadPart);
    return ticks / (double)t->frequency.QuadPart;
}

double timer_elapsed_ms(const Timer* t) {
    return timer_elapsed_s(t) * 1000.0;
}

double timer_elapsed_us(const Timer* t) {
    return timer_elapsed_s(t) * 1000000.0;
}
