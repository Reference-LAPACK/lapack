
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zsum_x, BLAS_ZSUM_X)
  (int *n, const void *x, int *incx, void *sum, int *prec) {
  BLAS_zsum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
