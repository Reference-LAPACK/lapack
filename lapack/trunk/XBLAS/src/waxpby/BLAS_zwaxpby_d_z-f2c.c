
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby_d_z(int n, const void *alpha, const double *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_zwaxpby_d_z, BLAS_ZWAXPBY_D_Z)
 
  (int *n, const void *alpha, const double *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_zwaxpby_d_z(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
