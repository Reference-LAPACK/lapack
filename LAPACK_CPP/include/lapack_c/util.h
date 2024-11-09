#ifndef LAPACK_C_UTIL_H
#define LAPACK_C_UTIL_H

// Define types for complex number
#ifdef __cplusplus
#include <complex>
#define lapack_float_complex std::complex<float>
#define lapack_double_complex std::complex<double>
#else
#define lapack_float_complex float _Complex
#define lapack_double_complex double _Complex
#endif // __cplusplus

// Default to 32 bit integer if not yet defined
#ifndef lapack_idx
#define lapack_idx int
#endif // lapack_idx

#endif