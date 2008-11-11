
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymv2_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const double *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);


extern void FC_FUNC_(blas_zsymv2_d_d_x, BLAS_ZSYMV2_D_D_X)
 
  (int *uplo, int *n, const void *alpha, const double *a, int *lda,
   const double *x_head, const double *x_tail, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_zsymv2_d_d_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		    *lda, x_head, x_tail, *incx, beta, y, *incy,
		    (enum blas_prec_type) *prec);
}
