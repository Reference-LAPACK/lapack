/*
 * c_dblas1.c
 *
 * The program is a C wrapper for dcblat1.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */
#include "cblas_test.h"
#include "cblas.h"
double F77_dasum(const CBLAS_INT *N, double *X, const CBLAS_INT *incX)
{
   return API_SUFFIX(cblas_dasum)(*N, X, *incX);
}

void F77_daxpy(const CBLAS_INT *N, const double *alpha, const double *X,
                    const CBLAS_INT *incX, double *Y, const CBLAS_INT *incY)
{
   API_SUFFIX(cblas_daxpy)(*N, *alpha, X, *incX, Y, *incY);
   return;
}

void F77_daxpby(const CBLAS_INT *N, const double *alpha, const double *X,
                    const CBLAS_INT *incX, const double *beta, double *Y, const CBLAS_INT *incY)
{
   API_SUFFIX(cblas_daxpby)(*N, *alpha, X, *incX, *beta, Y, *incY);
   return;
}


void F77_dcopy(const CBLAS_INT *N, double *X, const CBLAS_INT *incX,
                    double *Y, const CBLAS_INT *incY)
{
   API_SUFFIX(cblas_dcopy)(*N, X, *incX, Y, *incY);
   return;
}

double F77_ddot(const CBLAS_INT *N, const double *X, const CBLAS_INT *incX,
                const double *Y, const CBLAS_INT *incY)
{
   return API_SUFFIX(cblas_ddot)(*N, X, *incX, Y, *incY);
}

double F77_dnrm2(const CBLAS_INT *N, const double *X, const CBLAS_INT *incX)
{
   return API_SUFFIX(cblas_dnrm2)(*N, X, *incX);
}

void F77_drotg( double *a, double *b, double *c, double *s)
{
   API_SUFFIX(cblas_drotg)(a,b,c,s);
   return;
}

void F77_drot( const CBLAS_INT *N, double *X, const CBLAS_INT *incX, double *Y,
       const CBLAS_INT *incY, const double *c, const double *s)
{

   API_SUFFIX(cblas_drot)(*N,X,*incX,Y,*incY,*c,*s);
   return;
}

void F77_dscal(const CBLAS_INT *N, const double *alpha, double *X,
                         const CBLAS_INT *incX)
{
   API_SUFFIX(cblas_dscal)(*N, *alpha, X, *incX);
   return;
}

void F77_dswap( const CBLAS_INT *N, double *X, const CBLAS_INT *incX,
                          double *Y, const CBLAS_INT *incY)
{
   API_SUFFIX(cblas_dswap)(*N,X,*incX,Y,*incY);
   return;
}

CBLAS_INT F77_idamax(const CBLAS_INT *N, const double *X, const CBLAS_INT *incX)
{
   if (*N < 1 || *incX < 1) return(0);
   return (API_SUFFIX(cblas_idamax)(*N, X, *incX)+1);
}
