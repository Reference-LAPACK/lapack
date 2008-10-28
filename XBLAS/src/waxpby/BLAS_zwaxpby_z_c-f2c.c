
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby_z_c(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_zwaxpby_z_c, BLAS_ZWAXPBY_Z_C)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_zwaxpby_z_c(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
