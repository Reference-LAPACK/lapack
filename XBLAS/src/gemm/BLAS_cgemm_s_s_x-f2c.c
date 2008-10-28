
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemm_s_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const float *a, int lda,
		      const float *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cgemm_s_s_x, BLAS_CGEMM_S_S_X)
 
  (int *transa, int *transb, int *m, int *n, int *k, const void *alpha,
   const float *a, int *lda, const float *b, int *ldb, const void *beta,
   void *c, int *ldc, int *prec) {
  BLAS_cgemm_s_s_x(blas_colmajor, (enum blas_trans_type) *transa,
		   (enum blas_trans_type) *transb, *m, *n, *k, alpha, a, *lda,
		   b, *ldb, beta, c, *ldc, (enum blas_prec_type) *prec);
}
