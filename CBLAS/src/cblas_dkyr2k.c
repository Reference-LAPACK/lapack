/*
 *
 * cblas_dkyr2k.c
 * This program is a C interface to dkyr2k.
 * Written by Keita Teranishi
 * 4/6/1998
 *
 */

#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_dkyr2k)(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K,
                  const double alpha, const double  *A, const CBLAS_INT lda,
                  const double  *B, const CBLAS_INT ldb, const double beta,
                  double  *C, const CBLAS_INT ldc)
{
   char UL, TR;
   double minus_alpha;
#ifdef F77_CHAR
   F77_CHAR F77_TA, F77_UL;
#else
   #define F77_TR &TR
   #define F77_UL &UL
#endif

#ifdef F77_INT
   F77_INT F77_N=N, F77_K=K, F77_lda=lda, F77_ldb=ldb;
   F77_INT F77_ldc=ldc;
#else
   #define F77_N N
   #define F77_K K
   #define F77_lda lda
   #define F77_ldb ldb
   #define F77_ldc ldc
#endif

   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;
   CBLAS_CallFromC = 1;

   if( layout == CblasColMajor )
   {

      if( Uplo == CblasUpper) UL='U';
      else if ( Uplo == CblasLower ) UL='L';
      else
      {
         API_SUFFIX(cblas_xerbla)(2, "cblas_dkyr2k","Illegal Uplo setting, %d\n", Uplo);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Trans == CblasTrans) TR ='T';
      else if ( Trans == CblasConjTrans ) TR='C';
      else if ( Trans == CblasNoTrans )   TR='N';
      else
      {
         API_SUFFIX(cblas_xerbla)(3, "cblas_dkyr2k","Illegal Trans setting, %d\n", Trans);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }


      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_TR = C2F_CHAR(&TR);
      #endif

      F77_dkyr2k(F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda,
                      B, &F77_ldb, &beta, C, &F77_ldc);
   } else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if( Uplo == CblasUpper) UL='L';
      else if ( Uplo == CblasLower ) UL='U';
      else
      {
         API_SUFFIX(cblas_xerbla)(3, "cblas_dkyr2k","Illegal Uplo setting, %d\n", Uplo);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      if( Trans == CblasTrans) TR ='N';
      else if ( Trans == CblasConjTrans ) TR='N';
      else if ( Trans == CblasNoTrans )   TR='T';
      else
      {
         API_SUFFIX(cblas_xerbla)(3, "cblas_dkyr2k","Illegal Trans setting, %d\n", Trans);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_TR = C2F_CHAR(&TR);
      #endif

      minus_alpha = -alpha;
      F77_dkyr2k(F77_UL, F77_TR, &F77_N, &F77_K, &minus_alpha, A, &F77_lda, B,
                &F77_ldb, &beta, C, &F77_ldc);
   }
   else API_SUFFIX(cblas_xerbla)(1, "cblas_dkyr2k","Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
