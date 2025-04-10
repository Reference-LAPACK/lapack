/*
 *
 * cblas_skymv.c
 * This program is a C interface to skymv.
 * Written by Keita Teranishi
 * 4/6/1998
 *
 */

#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_skymv)(const CBLAS_LAYOUT layout,
                 const CBLAS_UPLO Uplo, const CBLAS_INT N,
                 const float alpha, const float  *A, const CBLAS_INT lda,
                 const float  *X, const CBLAS_INT incX, const float beta,
                 float  *Y, const CBLAS_INT incY)
{
   char UL;
   float minus_alpha;
#ifdef F77_CHAR
   F77_CHAR F77_UL;
#else
   #define F77_UL &UL
#endif
#ifdef F77_INT
   F77_INT F77_N=N, F77_lda=lda, F77_incX=incX, F77_incY=incY;
#else
   #define F77_N N
   #define F77_lda lda
   #define F77_incX incX
   #define F77_incY incY
#endif
   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      if (Uplo == CblasUpper) UL = 'U';
      else if (Uplo == CblasLower) UL = 'L';
      else
      {
         API_SUFFIX(cblas_xerbla)(2, "cblas_skymv","Illegal Uplo setting, %d\n",Uplo );
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
      #endif
      F77_skymv(F77_UL, &F77_N, &alpha, A, &F77_lda, X,
                     &F77_incX, &beta, Y, &F77_incY);
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if (Uplo == CblasUpper) UL = 'L';
      else if (Uplo == CblasLower) UL = 'U';
      else
      {
         API_SUFFIX(cblas_xerbla)(2, "cblas_skymv","Illegal Uplo setting, %d\n", Uplo);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
      #endif
      minus_alpha = -alpha;
      F77_skymv(F77_UL, &F77_N, &minus_alpha,
                     A ,&F77_lda, X,&F77_incX, &beta, Y, &F77_incY);
   }
   else API_SUFFIX(cblas_xerbla)(1, "cblas_skymv", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
