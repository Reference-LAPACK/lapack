
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgbmv2_s_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const float *a, int lda, const float *head_x,
		     const float *tail_x, int incx, const void *beta,
		     void *y, int incy);


extern void FC_FUNC_(blas_cgbmv2_s_s, BLAS_CGBMV2_S_S)
 
  (int *trans, int *m, int *n, int *kl, int *ku, const void *alpha,
   const float *a, int *lda, const float *head_x, const float *tail_x,
   int *incx, const void *beta, void *y, int *incy) {
  BLAS_cgbmv2_s_s(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl,
		  *ku, alpha, a, *lda, head_x, tail_x, *incx, beta, y, *incy);
}
