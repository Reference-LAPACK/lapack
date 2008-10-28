
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cwaxpby_s_s(int n, const void *alpha, const float *x, int incx,
		      const void *beta, const float *y, int incy, void *w,
		      int incw);


extern void FC_FUNC_(blas_cwaxpby_s_s, BLAS_CWAXPBY_S_S)
 
  (int *n, const void *alpha, const float *x, int *incx, const void *beta,
   const float *y, int *incy, void *w, int *incw) {
  BLAS_cwaxpby_s_s(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
