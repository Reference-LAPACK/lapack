
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zaxpby_d(int n, const void *alpha, const double *x, int incx,
		   const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_zaxpby_d, BLAS_ZAXPBY_D)
 
  (int *n, const void *alpha, const double *x, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_zaxpby_d(*n, alpha, x, *incx, beta, y, *incy);
}
