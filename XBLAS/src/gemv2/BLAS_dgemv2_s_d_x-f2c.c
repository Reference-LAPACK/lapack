
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemv2_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, double alpha, const float *a, int lda,
		       const double *head_x, const double *tail_x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);


extern void FC_FUNC_(blas_dgemv2_s_d_x, BLAS_DGEMV2_S_D_X)
 
  (int *trans, int *m, int *n, double *alpha, const float *a, int *lda,
   const double *head_x, const double *tail_x, int *incx, double *beta,
   double *y, int *incy, int *prec) {
  BLAS_dgemv2_s_d_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n,
		    *alpha, a, *lda, head_x, tail_x, *incx, *beta, y, *incy,
		    (enum blas_prec_type) *prec);
}
