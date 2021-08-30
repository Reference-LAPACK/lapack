/*
 * cblas_idamax.c
 *
 * The program is a C interface to idamax.
 * It calls the fortran wrapper before calling idamax.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
CBLAS_INDEX cblas_idamax( const CBLAS_INT N, const double *X, const CBLAS_INT incX)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_iamax;
#else
   #define F77_N N
   #define F77_incX incX
   CBLAS_INT F77_iamax;
#endif
   F77_idamax_sub( &F77_N, X, &F77_incX, &F77_iamax );
   return ( F77_iamax > 0 )
      ? (CBLAS_INDEX)( F77_iamax-1 )
      : (CBLAS_INDEX) 0;
}
