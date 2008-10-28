
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, const void *alpha,
		  const void *a, int lda, const void *x, int incx,
		  const void *beta, void *y, int incy,
		  enum blas_prec_type prec);


extern void FC_FUNC_(blas_cgbmv_x, BLAS_CGBMV_X)
 
  (int *trans, int *m, int *n, int *kl, int *ku, const void *alpha,
   const void *a, int *lda, const void *x, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_cgbmv_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl, *ku,
	       alpha, a, *lda, x, *incx, beta, y, *incy,
	       (enum blas_prec_type) *prec);
}
