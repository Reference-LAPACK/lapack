
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_daxpby_x(int n, double alpha, const double *x, int incx,
		   double beta, double *y,
		   int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_daxpby_x, BLAS_DAXPBY_X)
 
  (int *n, double *alpha, const double *x, int *incx, double *beta, double *y,
   int *incy, int *prec) {
  BLAS_daxpby_x(*n, *alpha, x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
