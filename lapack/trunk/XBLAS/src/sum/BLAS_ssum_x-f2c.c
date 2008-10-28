
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssum_x(int n, const float *x, int incx,
		 float *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ssum_x, BLAS_SSUM_X)
  (int *n, const float *x, int *incx, float *sum, int *prec) {
  BLAS_ssum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
