#ifndef MATGEN_H
#define MATGEN_H

#include <stddef.h>

double* matgen_random(size_t n, double low, double high, unsigned int seed);

double* matgen_random_vector(size_t n, double low, double high, unsigned int seed);

double* matgen_hilbert(size_t n);

double* matgen_rhs_from_exact(size_t n, const double* A, const double* x_exact);

double matgen_relative_error(size_t n, const double* x_approx, const double* x_exact);

#endif /* MATGEN_H */
