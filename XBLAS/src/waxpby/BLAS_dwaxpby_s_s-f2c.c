
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dwaxpby_s_s(int n, double alpha, const float *x, int incx,
		      double beta, const float *y, int incy, double *w,
		      int incw);


extern void FC_FUNC_(blas_dwaxpby_s_s, BLAS_DWAXPBY_S_S)
 
  (int *n, double *alpha, const float *x, int *incx, double *beta,
   const float *y, int *incy, double *w, int *incw) {
  BLAS_dwaxpby_s_s(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw);
}
