
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_caxpby_s(int n, const void *alpha, const float *x, int incx,
		   const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_caxpby_s, BLAS_CAXPBY_S)
 
  (int *n, const void *alpha, const float *x, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_caxpby_s(*n, alpha, x, *incx, beta, y, *incy);
}
