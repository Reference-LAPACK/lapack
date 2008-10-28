
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, const void *alpha, const void *a, int lda,
		   const void *head_x, const void *tail_x, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);


extern void FC_FUNC_(blas_cgemv2_x, BLAS_CGEMV2_X)
 
  (int *trans, int *m, int *n, const void *alpha, const void *a, int *lda,
   const void *head_x, const void *tail_x, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_cgemv2_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, alpha,
		a, *lda, head_x, tail_x, *incx, beta, y, *incy,
		(enum blas_prec_type) *prec);
}
