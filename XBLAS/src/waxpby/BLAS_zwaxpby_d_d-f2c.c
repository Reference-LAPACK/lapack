
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby_d_d(int n, const void *alpha, const double *x, int incx,
		      const void *beta, const double *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_zwaxpby_d_d, BLAS_ZWAXPBY_D_D)
 
  (int *n, const void *alpha, const double *x, int *incx, const void *beta,
   const double *y, int *incy, void *w, int *incw) {
  BLAS_zwaxpby_d_d(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
