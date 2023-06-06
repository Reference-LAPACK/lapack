/*
 * cblas_csrot.c
 *
 * The program is a C interface to csrot.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_csrot)(const CBLAS_INT N, void *X, const CBLAS_INT incX,
   void *Y, const CBLAS_INT incY, const float c, const float s)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
#endif
   F77_csrot(&F77_N, X, &F77_incX, Y, &F77_incY, &c, &s);
   return;
}
