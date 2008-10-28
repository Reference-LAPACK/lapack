
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, float alpha, const float *a, int lda,
		  const float *x, int incx, float beta, float *y,
		  int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_sgemv_x, BLAS_SGEMV_X)
 
  (int *trans, int *m, int *n, float *alpha, const float *a, int *lda,
   const float *x, int *incx, float *beta, float *y, int *incy, int *prec) {
  BLAS_sgemv_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *alpha,
	       a, *lda, x, *incx, *beta, y, *incy,
	       (enum blas_prec_type) *prec);
}
