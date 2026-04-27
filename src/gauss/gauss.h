#ifndef GAUSS_H
#define GAUSS_H

#include <stddef.h>
#include "../timer/timer.h"

#define ALGEBRA_OK        0
#define ALGEBRA_SINGULAR  1
#define ALGEBRA_ALLOC_ERR 2

#define GAUSS_METHOD_CLASSIC 0
#define GAUSS_METHOD_PIVOT   1

/* Hilbert matrices used in the lab can produce pivots below 1e-12 while still
   being nonzero in double precision. */
#define GAUSS_EPS 1e-16

double* gauss_alloc_matrix(size_t n);

void gauss_free_matrix(double* A);

double gauss_at(const double* A, size_t n, size_t i, size_t j);

double* gauss_at_ptr(double* A, size_t n, size_t i, size_t j);

void gauss_print_matrix(size_t n, const double* A, const char* label);

void gauss_print_vector(size_t n, const double* v, const char* label);

int gauss_solve(size_t n, const double* A, const double* b, double* x,
                int method, double* elapsed_ms);

int gauss_forward_substitution(size_t n, const double* L, const double* b,
                               double* y, double* elapsed_ms);

int gauss_back_substitution(size_t n, const double* U, const double* y,
                            double* x, double* elapsed_ms);

double gauss_residual(size_t n, const double* A, const double* b, const double* x);

#endif
