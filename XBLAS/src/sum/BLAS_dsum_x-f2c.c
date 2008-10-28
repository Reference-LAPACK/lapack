
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsum_x(int n, const double *x, int incx,
		 double *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dsum_x, BLAS_DSUM_X)
  (int *n, const double *x, int *incx, double *sum, int *prec) {
  BLAS_dsum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
