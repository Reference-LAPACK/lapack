/*
 * cblas_zgeru.c
 * The program is a C interface to zgeru.
 *
 * Keita Teranishi  5/20/98
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_zgeru)(const CBLAS_LAYOUT layout, const CBLAS_INT M, const CBLAS_INT N,
                 const void *alpha, const void *X, const CBLAS_INT incX,
                 const void *Y, const CBLAS_INT incY, void *A, const CBLAS_INT lda)
{
#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_lda=lda, F77_incX=incX, F77_incY=incY;
#else
   #define F77_M M
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
   #define F77_lda lda
#endif

   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;
   CBLAS_CallFromC = 1;

   if (layout == CblasColMajor)
   {
      F77_zgeru( &F77_M, &F77_N, alpha, X, &F77_incX, Y, &F77_incY, A,
                      &F77_lda);
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      F77_zgeru( &F77_N, &F77_M, alpha, Y, &F77_incY, X, &F77_incX, A,
                      &F77_lda);
   }
   else API_SUFFIX(cblas_xerbla)(1, "cblas_zgeru", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
