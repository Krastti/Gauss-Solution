#ifndef LU_H
#define LU_H

#include <stddef.h>
#include "../gauss/gauss.h"

int lu_decompose(size_t n, const double* A, double* L, double* U, double* elapsed_ms);

int lu_solve(size_t n, const double* L, const double* U, const double* b,
             double* x, double* elapsed_ms);

int lu_decompose_solve(size_t n, const double* A, const double* b, double* x,
                       double* decompose_ms, double* solve_ms, double* total_ms);

#endif
