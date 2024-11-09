#ifndef LAPACK_C_STEQR_HPP
#define LAPACK_C_STEQR_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

    void lapack_c_ssteqr(char compz,
                         lapack_idx n,
                         float *d,
                         float *e,
                         float *Z,
                         lapack_idx ldz,
                         float *work,
                         lapack_idx *info);

    void lapack_c_dsteqr(char compz,
                         lapack_idx n,
                         double *d,
                         double *e,
                         double *Z,
                         lapack_idx ldz,
                         double *work,
                         lapack_idx *info);

    void lapack_c_csteqr(char compz,
                         lapack_idx n,
                         float *d,
                         float *e,
                         lapack_float_complex *Z,
                         lapack_idx ldz,
                         float *work,
                         lapack_idx *info);

    void lapack_c_zsteqr(char compz,
                         lapack_idx n,
                         double *d,
                         double *e,
                         lapack_double_complex *Z,
                         lapack_idx ldz,
                         double *work,
                         lapack_idx *info);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPACK_C_STEQR_HPP