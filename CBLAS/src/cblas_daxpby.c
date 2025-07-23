/*
 * cblas_daxpby.c
 *
 * The program is a C interface to daxpby.
 *
 * Written by Martin Koehler.  08/26/2024
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_daxpby)( const CBLAS_INT N, const double alpha, const double *X,
                       const CBLAS_INT incX, const double beta, double *Y, const CBLAS_INT incY)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_daxpby( &F77_N, &alpha, X, &F77_incX, &beta, Y, &F77_incY);
}
