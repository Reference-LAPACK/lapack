
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dwaxpby_d_s(int n, double alpha, const double *x, int incx,
		      double beta, const float *y, int incy, double *w,
		      int incw);


extern void FC_FUNC_(blas_dwaxpby_d_s, BLAS_DWAXPBY_D_S)
 
  (int *n, double *alpha, const double *x, int *incx, double *beta,
   const float *y, int *incy, double *w, int *incw) {
  BLAS_dwaxpby_d_s(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw);
}
