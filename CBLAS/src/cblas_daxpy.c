/*
 * cblas_daxpy.c
 *
 * The program is a C interface to daxpy.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_daxpy)( const CBLAS_INT N, const double alpha, const double *X,
                       const CBLAS_INT incX, double *Y, const CBLAS_INT incY)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_daxpy( &F77_N, &alpha, X, &F77_incX, Y, &F77_incY);
}
