#include "gauss.h"
#include "../logging/logger.h"
#include "../timer/timer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double vector_at(const double* v, size_t n, size_t i) {
    if (v == NULL || n == 0 || i >= n) {
        LOG_ERROR("%s", "Некорректные параметры доступа к вектору!");
        return 0.0;
    }

    return v[i];
}

static double* vector_at_ptr(double* v, size_t n, size_t i) {
    if (v == NULL || n == 0 || i >= n) {
        LOG_ERROR("%s", "Некорректные параметры доступа к вектору!");
        return NULL;
    }

    return &v[i];
}

static int vector_set(double* v, size_t n, size_t i, double value) {
    double* element = vector_at_ptr(v, n, i);

    if (element == NULL) {
        return GAUSS_ERROR;
    }

    *element = value;
    return GAUSS_OK;
}

double gauss_at(const double* A, size_t n, size_t i, size_t j) {
    if (A == NULL || n == 0) {
        LOG_ERROR("%s", "Некорректные параметры доступа к матрице!");
        return 0.0;
    }

    return A[i * n + j];
}

double* gauss_at_ptr(double* A, size_t n, size_t i, size_t j) {
    if (A == NULL || n == 0) {
        LOG_ERROR("%s", "Некорректные параметры доступа к матрице!");
        return NULL;
    }

    return &A[i * n + j];
}

static int gauss_set(double* A, size_t n, size_t i, size_t j, double value) {
    double* element = gauss_at_ptr(A, n, i, j);

    if (element == NULL) {
        return GAUSS_ERROR;
    }

    *element = value;
    return GAUSS_OK;
}

double* gauss_alloc_matrix(size_t n) {
    if (n == 0) {
        LOG_ERROR("%s", "Размер матрицы должен быть больше нуля!");
        return NULL;
    }

    return calloc(n * n, sizeof(double));
}

void gauss_free_matrix(double* A) {
    free(A);
}

void gauss_print_vector(size_t n, const double* v, const char* label) {
    if (v == NULL) {
        LOG_ERROR("%s", "Попытка печати несуществующего вектора!");
        return;
    }
    
    printf("%s:\n| ", label);
    for (size_t i = 0; i < n; i++) {
        printf("%5.2f ", vector_at(v, n, i));
    }
    printf("|\n\n");
}

static int check_system(size_t n, const double* A, const double* b, const double* x) {
    if (n == 0) {
        LOG_ERROR("%s", "Размер системы должен быть больше нуля!");
        return 0;
    }

    if (A == NULL || b == NULL || x == NULL) {
        LOG_ERROR("%s", "Передан нулевой указатель в функцию решения СЛАУ!");
        return 0;
    }

    return 1;
}

static int swap_rows(double* A, double* b, size_t n, size_t row1, size_t row2) {
    double temp;

    for (size_t j = 0; j < n; j++) {
        temp = gauss_at(A, n, row1, j);
        if (gauss_set(A, n, row1, j, gauss_at(A, n, row2, j)) != GAUSS_OK ||
            gauss_set(A, n, row2, j, temp) != GAUSS_OK) {
            return GAUSS_ERROR;
        }
    }

    temp = vector_at(b, n, row1);
    if (vector_set(b, n, row1, vector_at(b, n, row2)) != GAUSS_OK ||
        vector_set(b, n, row2, temp) != GAUSS_OK) {
        return GAUSS_ERROR;
    }

    return GAUSS_OK;
}

static int forward_classic(size_t n, double* A, double* b) {
    double pivot;
    double factor;

    for (size_t col = 0; col < n; col++) {
        pivot = gauss_at(A, n, col, col);
        if (fabs(pivot) < GAUSS_EPS) {
            LOG_ERROR("Нулевой диагональный элемент на шаге %zu", col);
            return GAUSS_ERROR;
        }

        for (size_t row = col + 1; row < n; row++) {
            factor = gauss_at(A, n, row, col) / pivot;
            if (vector_set(b, n, row, vector_at(b, n, row) - factor * vector_at(b, n, col)) != GAUSS_OK) {
                return GAUSS_ERROR;
            }

            for (size_t j = col; j < n; j++) {
                if (gauss_set(A, n, row, j, gauss_at(A, n, row, j) - factor * gauss_at(A, n, col, j)) != GAUSS_OK) {
                    return GAUSS_ERROR;
                }
            }
            if (gauss_set(A, n, row, col, 0.0) != GAUSS_OK) {
                return GAUSS_ERROR;
            }
        }
    }

    return GAUSS_OK;
}

static int forward_pivot(size_t n, double* A, double* b) {
    double pivot;
    double factor;

    for (size_t col = 0; col < n; col++) {
        size_t pivot_row = col;
        double max_value = fabs(gauss_at(A, n, col, col));

        for (size_t row = col + 1; row < n; row++) {
            double value = fabs(gauss_at(A, n, row, col));
            if (value > max_value) {
                max_value = value;
                pivot_row = row;
            }
        }

        if (max_value < GAUSS_EPS) {
            LOG_ERROR("Ведущий элемент на шаге %zu слишком мал", col);
            return GAUSS_ERROR;
        }

        if (pivot_row != col) {
            if (swap_rows(A, b, n, col, pivot_row) != GAUSS_OK) {
                return GAUSS_ERROR;
            }
        }

        pivot = gauss_at(A, n, col, col);
        for (size_t row = col + 1; row < n; row++) {
            factor = gauss_at(A, n, row, col) / pivot;
            if (vector_set(b, n, row, vector_at(b, n, row) - factor * vector_at(b, n, col)) != GAUSS_OK) {
                return GAUSS_ERROR;
            }

            for (size_t j = col; j < n; j++) {
                if (gauss_set(A, n, row, j, gauss_at(A, n, row, j) - factor * gauss_at(A, n, col, j)) != GAUSS_OK) {
                    return GAUSS_ERROR;
                }
            }
            if (gauss_set(A, n, row, col, 0.0) != GAUSS_OK) {
                return GAUSS_ERROR;
            }
        }
    }

    return GAUSS_OK;
}

static int back_substitution(size_t n, const double* A, const double* b, double* x) {
    size_t i = n;

    while (i > 0) {
        i--;

        if (fabs(gauss_at(A, n, i, i)) < GAUSS_EPS) {
            LOG_ERROR("Нулевой диагональный элемент при обратной подстановке: строка %zu", i);
            return GAUSS_ERROR;
        }

        if (vector_set(x, n, i, vector_at(b, n, i)) != GAUSS_OK) {
            return GAUSS_ERROR;
        }
        for (size_t j = i + 1; j < n; j++) {
            if (vector_set(x, n, i, vector_at(x, n, i) - gauss_at(A, n, i, j) * vector_at(x, n, j)) != GAUSS_OK) {
                return GAUSS_ERROR;
            }
        }
        if (vector_set(x, n, i, vector_at(x, n, i) / gauss_at(A, n, i, i)) != GAUSS_OK) {
            return GAUSS_ERROR;
        }
    }

    return GAUSS_OK;
}

int gauss_solve(size_t n, const double* A, const double* b, double* x, int method, double* elapsed_ms) {
    Timer timer;
    double* work_A;
    double* work_b;
    int status;

    if (!check_system(n, A, b, x)) {
        return GAUSS_ERROR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    work_A = malloc(n * n * sizeof(double));
    work_b = malloc(n * sizeof(double));
    if (work_A == NULL || work_b == NULL) {
        free(work_A);
        free(work_b);
        return GAUSS_ERROR;
    }

    memcpy(work_A, A, n * n * sizeof(double));
    memcpy(work_b, b, n * sizeof(double));

    if (method == GAUSS_METHOD_CLASSIC) {
        status = forward_classic(n, work_A, work_b);
    } else {
        status = forward_pivot(n, work_A, work_b);
    }

    if (status == GAUSS_OK) {
        status = back_substitution(n, work_A, work_b, x);
    }

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    free(work_A);
    free(work_b);
    return status;
}

int gauss_forward_substitution(size_t n, const double* L, const double* b, double* y, double* elapsed_ms) {
    Timer timer;

    if (!check_system(n, L, b, y)) {
        return GAUSS_ERROR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    for (size_t i = 0; i < n; i++) {
        if (fabs(gauss_at(L, n, i, i)) < GAUSS_EPS) {
            if (elapsed_ms != NULL) {
                timer_stop(&timer);
                *elapsed_ms = timer_elapsed_ms(&timer);
            }
            return GAUSS_ERROR;
        }

        if (vector_set(y, n, i, vector_at(b, n, i)) != GAUSS_OK) {
            return GAUSS_ERROR;
        }
        for (size_t j = 0; j < i; j++) {
            if (vector_set(y, n, i, vector_at(y, n, i) - gauss_at(L, n, i, j) * vector_at(y, n, j)) != GAUSS_OK) {
                return GAUSS_ERROR;
            }
        }
        if (vector_set(y, n, i, vector_at(y, n, i) / gauss_at(L, n, i, i)) != GAUSS_OK) {
            return GAUSS_ERROR;
        }
    }

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    return GAUSS_OK;
}

int gauss_back_substitution(size_t n, const double* U, const double* y, double* x, double* elapsed_ms) {
    Timer timer;
    int status;

    if (!check_system(n, U, y, x)) {
        return GAUSS_ERROR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    status = back_substitution(n, U, y, x);

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    return status;
}

double gauss_residual(size_t n, const double* A, const double* b, const double* x) {
    double sum = 0.0;

    if (!check_system(n, A, b, x)) {
        return -1.0;
    }

    for (size_t i = 0; i < n; i++) {
        double ax = 0.0;

        for (size_t j = 0; j < n; j++) {
            ax = ax + gauss_at(A, n, i, j) * vector_at(x, n, j);
        }

        sum = sum + (ax - vector_at(b, n, i)) * (ax - vector_at(b, n, i));
    }

    return sqrt(sum);
}
