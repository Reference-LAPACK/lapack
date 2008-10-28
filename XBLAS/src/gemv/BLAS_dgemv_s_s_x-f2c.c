
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemv_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double alpha, const float *a, int lda,
		      const float *x, int incx, double beta, double *y,
		      int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dgemv_s_s_x, BLAS_DGEMV_S_S_X)
 
  (int *trans, int *m, int *n, double *alpha, const float *a, int *lda,
   const float *x, int *incx, double *beta, double *y, int *incy, int *prec) {
  BLAS_dgemv_s_s_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n,
		   *alpha, a, *lda, x, *incx, *beta, y, *incy,
		   (enum blas_prec_type) *prec);
}
