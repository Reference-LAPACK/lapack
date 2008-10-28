
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_daxpby_s(int n, double alpha, const float *x, int incx,
		   double beta, double *y, int incy);


extern void FC_FUNC_(blas_daxpby_s, BLAS_DAXPBY_S)
 
  (int *n, double *alpha, const float *x, int *incx, double *beta, double *y,
   int *incy) {
  BLAS_daxpby_s(*n, *alpha, x, *incx, *beta, y, *incy);
}
