
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_daxpby_s_x(int n, double alpha, const float *x, int incx,
		     double beta, double *y,
		     int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_daxpby_s_x, BLAS_DAXPBY_S_X)
 
  (int *n, double *alpha, const float *x, int *incx, double *beta, double *y,
   int *incy, int *prec) {
  BLAS_daxpby_s_x(*n, *alpha, x, *incx, *beta, y, *incy,
		  (enum blas_prec_type) *prec);
}
