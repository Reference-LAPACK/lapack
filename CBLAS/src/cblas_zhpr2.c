/*
 * cblas_zhpr2.c
 * The program is a C interface to zhpr2.
 *
 * Keita Teranishi  5/20/98
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_zhpr2)(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                      const CBLAS_INT N,const void *alpha, const void *X,
                      const CBLAS_INT incX,const void *Y, const CBLAS_INT incY, void *Ap)

{
   char UL;
#ifdef F77_CHAR
   F77_CHAR F77_UL;
#else
   #define F77_UL &UL
#endif

#ifdef F77_INT
   F77_INT F77_N=N,  F77_incX=incX, F77_incY=incY;
#else
   CBLAS_INT incx = incX, incy = incY;
   #define F77_N N
   #define F77_incX incx
   #define F77_incY incy
#endif
   CBLAS_INT n, i, j;
   double *x, *xx, *y,
         *yy, *stx, *sty;

   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   memcpy(&x,&X,sizeof(double*));
   memcpy(&xx,&X,sizeof(double*));
   memcpy(&y,&Y,sizeof(double*));
   memcpy(&yy,&Y,sizeof(double*));

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      if (Uplo == CblasLower) UL = 'L';
      else if (Uplo == CblasUpper) UL = 'U';
      else
      {
         API_SUFFIX(cblas_xerbla)(2, "cblas_zhpr2","Illegal Uplo setting, %d\n",Uplo );
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
      #endif

      F77_zhpr2(F77_UL, &F77_N, alpha, X, &F77_incX, Y, &F77_incY, Ap);

   }  else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if (Uplo == CblasUpper) UL = 'L';
      else if (Uplo == CblasLower) UL = 'U';
      else
      {
         API_SUFFIX(cblas_xerbla)(2, "cblas_zhpr2","Illegal Uplo setting, %d\n", Uplo);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
      #endif
      if (N > 0)
      {
         n = N << 1;
         x = malloc(n*sizeof(double));
         y = malloc(n*sizeof(double));
         stx = x + n;
         sty = y + n;
         if( incX > 0 )
            i = incX << 1;
         else
            i = incX *(-2);

         if( incY > 0 )
            j = incY << 1;
         else
            j = incY *(-2);
         do
         {
            *x = *xx;
            x[1] = -xx[1];
            x += 2;
            xx += i;
         } while (x != stx);
         do
         {
            *y = *yy;
            y[1] = -yy[1];
            y += 2;
            yy += j;
         }
         while (y != sty);
         x -= n;
         y -= n;

         #ifdef F77_INT
            if(incX > 0 )
               F77_incX = 1;
            else
               F77_incX = -1;

            if(incY > 0 )
               F77_incY = 1;
            else
               F77_incY = -1;

         #else
            if(incX > 0 )
               incx = 1;
            else
               incx = -1;

            if(incY > 0 )
               incy = 1;
            else
               incy = -1;
         #endif

      }  else
      {

        memcpy(&x,&X,sizeof(double*));
        memcpy(&y,&Y,sizeof(double*));
      }
      F77_zhpr2(F77_UL, &F77_N, alpha, y, &F77_incY, x, &F77_incX, Ap);
   }
   else
   {
      API_SUFFIX(cblas_xerbla)(1, "cblas_zhpr2","Illegal layout setting, %d\n", layout);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
   }
   if(X!=x)
      free(x);
   if(Y!=y)
      free(y);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
