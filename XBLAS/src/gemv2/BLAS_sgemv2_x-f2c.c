
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, float alpha, const float *a, int lda,
		   const float *head_x, const float *tail_x, int incx,
		   float beta, float *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_sgemv2_x, BLAS_SGEMV2_X)
 
  (int *trans, int *m, int *n, float *alpha, const float *a, int *lda,
   const float *head_x, const float *tail_x, int *incx, float *beta, float *y,
   int *incy, int *prec) {
  BLAS_sgemv2_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *alpha,
		a, *lda, head_x, tail_x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
