#include "matgen.h"
#include "../gauss/gauss.h"
#include "../logging/logger.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

static double random_0_1(unsigned int* seed) {
    *seed = 1664525u * (*seed) + 1013904223u;
    return (double)(*seed) / (double)UINT_MAX;
}

double* matgen_random(size_t n, double low, double high, unsigned int seed) {
    double* A;
    double value;

    if (n == 0 || low >= high) {
        LOG_ERROR("%s", "Некорректные параметры matgen_random!");
        return NULL;
    }

    A = gauss_alloc_matrix(n);
    if (A == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            value = low + (high - low) * random_0_1(&seed);
            A[i * n + j] = value;
        }
    }

    return A;
}

double* matgen_random_vector(size_t n, double low, double high, unsigned int seed) {
    double* v;

    if (n == 0 || low >= high) {
        LOG_ERROR("%s", "Некорректные параметры matgen_random_vector!");
        return NULL;
    }

    v = malloc(n * sizeof(double));
    if (v == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        v[i] = low + (high - low) * random_0_1(&seed);
    }

    return v;
}

double* matgen_hilbert(size_t n) {
    double* H;

    if (n == 0) {
        LOG_ERROR("%s", "Размер матрицы Гильберта должен быть больше нуля!");
        return NULL;
    }

    H = gauss_alloc_matrix(n);
    if (H == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            H[i * n + j] = 1.0 / (double)(i + j + 1);
        }
    }

    return H;
}

double* matgen_rhs_from_exact(size_t n, const double* A, const double* x_exact) {
    double* b;

    if (n == 0 || A == NULL || x_exact == NULL) {
        LOG_ERROR("%s", "Некорректные параметры matgen_rhs_from_exact!");
        return NULL;
    }

    b = malloc(n * sizeof(double));
    if (b == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        b[i] = 0.0;
        for (size_t j = 0; j < n; j++) {
            b[i] = b[i] + A[i * n + j] * x_exact[j];
        }
    }

    return b;
}

double matgen_relative_error(size_t n, const double* x_approx, const double* x_exact) {
    double diff_norm = 0.0;
    double exact_norm = 0.0;

    if (n == 0 || x_approx == NULL || x_exact == NULL) {
        LOG_ERROR("%s", "Некорректные параметры matgen_relative_error!");
        return -1.0;
    }

    for (size_t i = 0; i < n; i++) {
        double diff = x_approx[i] - x_exact[i];
        diff_norm = diff_norm + diff * diff;
        exact_norm = exact_norm + x_exact[i] * x_exact[i];
    }

    if (exact_norm == 0.0) {
        return -1.0;
    }

    return sqrt(diff_norm) / sqrt(exact_norm);
}
