#ifndef LAPACK_C_ROT_HPP
#define LAPACK_C_ROT_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

    void lapack_c_ssteqr3(char compz,
                         lapack_idx n,
                         float *d,
                         float *e,
                         float *Z,
                         lapack_idx ldz,
                         float *work,
                         lapack_idx lwork,
                         float *rwork,
                         lapack_idx lrwork,
                         lapack_idx *info);

    void lapack_c_dsteqr3(char compz,
                         lapack_idx n,
                         double *d,
                         double *e,
                         double *Z,
                         lapack_idx ldz,
                         double *work,
                         lapack_idx lwork,
                         double *rwork,
                         lapack_idx lrwork,
                         lapack_idx *info);

    void lapack_c_csteqr3(char compz,
                         lapack_idx n,
                         float *d,
                         float *e,
                         lapack_float_complex *Z,
                         lapack_idx ldz,
                         lapack_float_complex *work,
                         lapack_idx lwork,
                         float *rwork,
                         lapack_idx lrwork,
                         lapack_idx *info);

    void lapack_c_zsteqr3(char compz,
                         lapack_idx n,
                         double *d,
                         double *e,
                         lapack_double_complex *Z,
                         lapack_idx ldz,
                         lapack_double_complex *work,
                         lapack_idx lwork,
                         double *rwork,
                         lapack_idx lrwork,
                         lapack_idx *info);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif