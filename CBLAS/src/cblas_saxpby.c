/*
 * cblas_saxpby.c
 *
 * The program is a C interface to saxpby.
 * It calls the fortran wrapper before calling saxpby.
 *
 * Written by Martin Koehler, 08/24/2024
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_saxpby)( const CBLAS_INT N, const float alpha, const float *X,
                       const CBLAS_INT incX, const float beta, float *Y, const CBLAS_INT incY)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_saxpby( &F77_N, &alpha, X, &F77_incX, &beta, Y, &F77_incY);
}
