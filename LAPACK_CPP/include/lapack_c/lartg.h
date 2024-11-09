#ifndef LAPACK_C_LARTG_HPP
#define LAPACK_C_LARTG_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_slartg(float f, float g, float* cs, float* sn, float* r);

void lapack_c_dlartg(double f, double g, double* cs, double* sn, double* r);

void lapack_c_clartg(lapack_float_complex f, lapack_float_complex g, float* cs, lapack_float_complex* sn, lapack_float_complex* r);

void lapack_c_zlartg(lapack_double_complex f, lapack_double_complex g, double* cs, lapack_double_complex* sn, lapack_double_complex* r);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif