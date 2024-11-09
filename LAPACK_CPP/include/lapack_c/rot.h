#ifndef LAPACK_C_ROT_HPP
#define LAPACK_C_ROT_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_drot(lapack_idx n,
                    double* x,
                    lapack_idx incx,
                    double* y,
                    lapack_idx incy,
                    double c,
                    double s);

void lapack_c_srot(lapack_idx n,
                    float* x,
                    lapack_idx incx,
                    float* y,
                    lapack_idx incy,
                    float c,
                    float s);

void lapack_c_crot(lapack_idx n,
                    lapack_float_complex* x,
                    lapack_idx incx,
                    lapack_float_complex* y,
                    lapack_idx incy,
                    float c,
                    lapack_float_complex s);

void lapack_c_zrot(lapack_idx n,
                    lapack_double_complex* x,
                    lapack_idx incx,
                    lapack_double_complex* y,
                    lapack_idx incy,
                    double c,
                    lapack_double_complex s);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif