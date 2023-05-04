/*
 * cblas_zdrot.c
 *
 * The program is a C interface to zdrot.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_zdrot)(const CBLAS_INT N, void *X, const CBLAS_INT incX,
   void *Y, const CBLAS_INT incY, const double c, const double s)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_zdrot(&F77_N, X, &F77_incX, Y, &F77_incY, &c, &s);
   return;
}
