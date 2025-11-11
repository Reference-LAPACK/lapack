/*
 * cblas_caxpby.c
 *
 * The program is a C interface to caxpby.
 *
 * Written by Martin Koehler.  08/26/2024
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_caxpby)( const CBLAS_INT N, const void *alpha, const void *X,
                       const CBLAS_INT incX, const void *beta, void *Y, const CBLAS_INT incY)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_caxpby( &F77_N, alpha, X, &F77_incX, beta, Y, &F77_incY);
}
