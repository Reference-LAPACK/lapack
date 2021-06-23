/*
 * cblas_izamax.c
 *
 * The program is a C interface to izamax.
 * It calls the fortran wrapper before calling izamax.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
CBLAS_INDEX cblas_izamax( const CBLAS_INT N, const void *X, const CBLAS_INT incX)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_iamax;
#else
   #define F77_N N
   #define F77_incX incX
   CBLAS_INT F77_iamax;
#endif
   F77_izamax_sub( &F77_N, X, &F77_incX, &F77_iamax );
   return ( F77_iamax > 0 )
      ? (CBLAS_INDEX)( F77_iamax-1 )
      : (CBLAS_INDEX) 0;
}
