
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgbmv_s_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, double alpha,
		    const float *a, int lda, const float *x, int incx,
		    double beta, double *y, int incy);


extern void FC_FUNC_(blas_dgbmv_s_s, BLAS_DGBMV_S_S)
 
  (int *trans, int *m, int *n, int *kl, int *ku, double *alpha,
   const float *a, int *lda, const float *x, int *incx, double *beta,
   double *y, int *incy) {
  BLAS_dgbmv_s_s(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl,
		 *ku, *alpha, a, *lda, x, *incx, *beta, y, *incy);
}
