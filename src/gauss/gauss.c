#include "gauss.h"
#include "../logging/logger.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void gauss_print_matrix(size_t n, const double* A, const char* label) {
    if (A == NULL) {
        LOG_ERROR("%s", "Попытка печати несуществующей матрицы!");
        return;
    }

    if (label == NULL) {
        label = "Матрица";
    }

    printf("%s %zux%zu:\n", label, n, n);
    for (size_t i = 0; i < n; i++) {
        printf("| ");
        for (size_t j = 0; j < n; j++) {
            printf("%10.4f ", A[i * n + j]);
        }
        printf("|\n");
    }
    printf("\n");
}

void gauss_print_vector(size_t n, const double* v, const char* label) {
    if (v == NULL) {
        LOG_ERROR("%s", "Попытка печати несуществующего вектора!");
        return;
    }

    if (label == NULL) {
        label = "Вектор";
    }

    printf("%s:\n| ", label);
    for (size_t i = 0; i < n; i++) {
        printf("%10.4f ", v[i]);
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

static void swap_rows(double* A, double* b, size_t n, size_t row1, size_t row2) {
    double temp;

    for (size_t j = 0; j < n; j++) {
        temp = A[row1 * n + j];
        A[row1 * n + j] = A[row2 * n + j];
        A[row2 * n + j] = temp;
    }

    temp = b[row1];
    b[row1] = b[row2];
    b[row2] = temp;
}

static int forward_classic(size_t n, double* A, double* b) {
    double pivot;
    double factor;

    for (size_t col = 0; col < n; col++) {
        pivot = A[col * n + col];
        if (fabs(pivot) < GAUSS_EPS) {
            LOG_ERROR("Нулевой диагональный элемент на шаге %zu", col);
            return ALGEBRA_SINGULAR;
        }

        for (size_t row = col + 1; row < n; row++) {
            factor = A[row * n + col] / pivot;
            b[row] = b[row] - factor * b[col];

            for (size_t j = col; j < n; j++) {
                A[row * n + j] = A[row * n + j] - factor * A[col * n + j];
            }
            A[row * n + col] = 0.0;
        }
    }

    return ALGEBRA_OK;
}

static int forward_pivot(size_t n, double* A, double* b) {
    double pivot;
    double factor;

    for (size_t col = 0; col < n; col++) {
        size_t pivot_row = col;
        double max_value = fabs(A[col * n + col]);

        for (size_t row = col + 1; row < n; row++) {
            double value = fabs(A[row * n + col]);
            if (value > max_value) {
                max_value = value;
                pivot_row = row;
            }
        }

        if (max_value < GAUSS_EPS) {
            LOG_ERROR("Ведущий элемент на шаге %zu слишком мал", col);
            return ALGEBRA_SINGULAR;
        }

        if (pivot_row != col) {
            swap_rows(A, b, n, col, pivot_row);
        }

        pivot = A[col * n + col];
        for (size_t row = col + 1; row < n; row++) {
            factor = A[row * n + col] / pivot;
            b[row] = b[row] - factor * b[col];

            for (size_t j = col; j < n; j++) {
                A[row * n + j] = A[row * n + j] - factor * A[col * n + j];
            }
            A[row * n + col] = 0.0;
        }
    }

    return ALGEBRA_OK;
}

static int back_substitution(size_t n, const double* A, const double* b, double* x) {
    size_t i = n;

    while (i > 0) {
        i--;

        if (fabs(A[i * n + i]) < GAUSS_EPS) {
            LOG_ERROR("Нулевой диагональный элемент при обратной подстановке: строка %zu", i);
            return ALGEBRA_SINGULAR;
        }

        x[i] = b[i];
        for (size_t j = i + 1; j < n; j++) {
            x[i] = x[i] - A[i * n + j] * x[j];
        }
        x[i] = x[i] / A[i * n + i];
    }

    return ALGEBRA_OK;
}

int gauss_solve(size_t n, const double* A, const double* b, double* x, int method, double* elapsed_ms) {
    Timer timer;
    double* work_A;
    double* work_b;
    int status;

    if (!check_system(n, A, b, x)) {
        return ALGEBRA_ALLOC_ERR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    work_A = malloc(n * n * sizeof(double));
    work_b = malloc(n * sizeof(double));
    if (work_A == NULL || work_b == NULL) {
        free(work_A);
        free(work_b);
        return ALGEBRA_ALLOC_ERR;
    }

    memcpy(work_A, A, n * n * sizeof(double));
    memcpy(work_b, b, n * sizeof(double));

    if (method == GAUSS_METHOD_CLASSIC) {
        status = forward_classic(n, work_A, work_b);
    } else {
        status = forward_pivot(n, work_A, work_b);
    }

    if (status == ALGEBRA_OK) {
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
        return ALGEBRA_ALLOC_ERR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    for (size_t i = 0; i < n; i++) {
        if (fabs(L[i * n + i]) < GAUSS_EPS) {
            if (elapsed_ms != NULL) {
                timer_stop(&timer);
                *elapsed_ms = timer_elapsed_ms(&timer);
            }
            return ALGEBRA_SINGULAR;
        }

        y[i] = b[i];
        for (size_t j = 0; j < i; j++) {
            y[i] = y[i] - L[i * n + j] * y[j];
        }
        y[i] = y[i] / L[i * n + i];
    }

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    return ALGEBRA_OK;
}

int gauss_back_substitution(size_t n, const double* U, const double* y, double* x, double* elapsed_ms) {
    Timer timer;
    int status;

    if (!check_system(n, U, y, x)) {
        return ALGEBRA_ALLOC_ERR;
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
            ax = ax + A[i * n + j] * x[j];
        }

        sum = sum + (ax - b[i]) * (ax - b[i]);
    }

    return sqrt(sum);
}
