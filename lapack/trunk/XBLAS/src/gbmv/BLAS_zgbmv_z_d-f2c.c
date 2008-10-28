
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zgbmv_z_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const double *x, int incx,
		    const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_zgbmv_z_d, BLAS_ZGBMV_Z_D)
 
  (int *trans, int *m, int *n, int *kl, int *ku, const void *alpha,
   const void *a, int *lda, const double *x, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_zgbmv_z_d(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl,
		 *ku, alpha, a, *lda, x, *incx, beta, y, *incy);
}
