
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cwaxpby_s_c(int n, const void *alpha, const float *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_cwaxpby_s_c, BLAS_CWAXPBY_S_C)
 
  (int *n, const void *alpha, const float *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_cwaxpby_s_c(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
