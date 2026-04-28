#include "lu.h"
#include "../gauss/gauss.h"
#include "../logging/logger.h"
#include "../timer/timer.h"

#include <math.h>
#include <stdlib.h>

static int check_lu_args(size_t n, const double* A, const char* name) {
    if (n == 0) {
        LOG_ERROR("%s: размер системы должен быть больше нуля!", name);
        return 0;
    }

    if (A == NULL) {
        LOG_ERROR("%s: передан нулевой указатель!", name);
        return 0;
    }

    return 1;
}

int lu_decompose(size_t n, const double* A, double* L, double* U, double* elapsed_ms) {
    Timer timer;

    if (!check_lu_args(n, A, "lu_decompose") || L == NULL || U == NULL) {
        return LU_ERROR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    for (size_t i = 0; i < n * n; i++) {
        L[i] = 0.0;
        U[i] = 0.0;
    }

    for (size_t i = 0; i < n; i++) {
        L[i * n + i] = 1.0;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t k = i; k < n; k++) {
            double sum = 0.0;

            for (size_t j = 0; j < i; j++) {
                sum = sum + L[i * n + j] * U[j * n + k];
            }

            U[i * n + k] = A[i * n + k] - sum;
        }

        if (fabs(U[i * n + i]) < LU_EPS) {
            if (elapsed_ms != NULL) {
                timer_stop(&timer);
                *elapsed_ms = timer_elapsed_ms(&timer);
            }
            return LU_ERROR;
        }

        for (size_t k = i + 1; k < n; k++) {
            double sum = 0.0;

            for (size_t j = 0; j < i; j++) {
                sum = sum + L[k * n + j] * U[j * n + i];
            }

            L[k * n + i] = (A[k * n + i] - sum) / U[i * n + i];
        }
    }

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    return LU_OK;
}

int lu_solve(size_t n, const double* L, const double* U, const double* b, double* x, double* elapsed_ms) {
    Timer timer;
    double* y;
    int status;

    if (n == 0 || L == NULL || U == NULL || b == NULL || x == NULL) {
        LOG_ERROR("%s", "Некорректные параметры lu_solve!");
        return LU_ERROR;
    }

    if (elapsed_ms != NULL) {
        timer_start(&timer);
    }

    y = malloc(n * sizeof(double));
    if (y == NULL) {
        return LU_ERROR;
    }

    status = gauss_forward_substitution(n, L, b, y, NULL);
    if (status == GAUSS_OK) {
        status = gauss_back_substitution(n, U, y, x, NULL);
    }

    if (elapsed_ms != NULL) {
        timer_stop(&timer);
        *elapsed_ms = timer_elapsed_ms(&timer);
    }

    free(y);
    return status;
}

int lu_decompose_solve(size_t n, const double* A, const double* b, double* x, double* decompose_ms, double* solve_ms, double* total_ms) {
    Timer timer;
    double* L;
    double* U;
    int status;

    if (n == 0 || A == NULL || b == NULL || x == NULL) {
        LOG_ERROR("%s", "Некорректные параметры lu_decompose_solve!");
        return LU_ERROR;
    }

    if (total_ms != NULL) {
        timer_start(&timer);
    }

    L = gauss_alloc_matrix(n);
    U = gauss_alloc_matrix(n);
    if (L == NULL || U == NULL) {
        gauss_free_matrix(L);
        gauss_free_matrix(U);
        return LU_ERROR;
    }

    status = lu_decompose(n, A, L, U, decompose_ms);
    if (status == LU_OK) {
        status = lu_solve(n, L, U, b, x, solve_ms);
    }

    if (total_ms != NULL) {
        timer_stop(&timer);
        *total_ms = timer_elapsed_ms(&timer);
    }

    gauss_free_matrix(L);
    gauss_free_matrix(U);
    return status;
}
