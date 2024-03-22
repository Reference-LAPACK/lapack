#ifndef LAPACK_C_LAEV2_HPP
#define LAPACK_C_LAEV2_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

    void lapack_c_slaev2(float a, float b, float c, float *rt1, float *rt2, float *cs1, float *sn1);

    void lapack_c_dlaev2(double a, double b, double c, double *rt1, double *rt2, double *cs1, double *sn1);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPACK_C_LAEV2_HPP