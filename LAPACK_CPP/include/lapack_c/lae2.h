#ifndef LAPACK_C_LAE2_HPP
#define LAPACK_C_LAE2_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

    void lapack_c_slae2(float a, float b, float c, float *rt1, float *rt2);

    void lapack_c_dlae2(double a, double b, double c, double *rt1, double *rt2);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPACK_C_LAE2_HPP