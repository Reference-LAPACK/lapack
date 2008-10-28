
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby_c_z(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_zwaxpby_c_z, BLAS_ZWAXPBY_C_Z)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_zwaxpby_c_z(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
