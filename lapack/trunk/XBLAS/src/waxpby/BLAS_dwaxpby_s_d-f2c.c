
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dwaxpby_s_d(int n, double alpha, const float *x, int incx,
		      double beta, const double *y, int incy, double *w,
		      int incw);


extern void FC_FUNC_(blas_dwaxpby_s_d, BLAS_DWAXPBY_S_D)
 
  (int *n, double *alpha, const float *x, int *incx, double *beta,
   const double *y, int *incy, double *w, int *incw) {
  BLAS_dwaxpby_s_d(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw);
}
