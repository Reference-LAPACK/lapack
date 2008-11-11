
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv2_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double alpha, const float *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);


extern void FC_FUNC_(blas_dsymv2_s_s_x, BLAS_DSYMV2_S_S_X)
 
  (int *uplo, int *n, double *alpha, const float *a, int *lda,
   const float *x_head, const float *x_tail, int *incx, double *beta,
   double *y, int *incy, int *prec) {
  BLAS_dsymv2_s_s_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		    *lda, x_head, x_tail, *incx, *beta, y, *incy,
		    (enum blas_prec_type) *prec);
}
