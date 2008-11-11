
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);


extern void FC_FUNC_(blas_zsymv2_x, BLAS_ZSYMV2_X)
 
  (int *uplo, int *n, const void *alpha, const void *a, int *lda,
   const void *x_head, const void *x_tail, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_zsymv2_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		*lda, x_head, x_tail, *incx, beta, y, *incy,
		(enum blas_prec_type) *prec);
}
