#ifndef LAPACK_C_LASRT_HPP
#define LAPACK_C_LASRT_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

    void lapack_c_slasrt(char id, lapack_idx n, float* d, lapack_idx *info);

    void lapack_c_dlasrt(char id, lapack_idx n, double* d, lapack_idx *info);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPACK_C_LASRT_HPP