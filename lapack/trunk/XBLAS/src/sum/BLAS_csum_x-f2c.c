
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_csum_x, BLAS_CSUM_X)
  (int *n, const void *x, int *incx, void *sum, int *prec) {
  BLAS_csum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
