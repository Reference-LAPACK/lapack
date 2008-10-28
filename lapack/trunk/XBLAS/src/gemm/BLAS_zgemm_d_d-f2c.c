
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zgemm_d_d(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const double *a, int lda,
		    const double *b, int ldb, const void *beta, void *c,
		    int ldc);


extern void FC_FUNC_(blas_zgemm_d_d, BLAS_ZGEMM_D_D)
 
  (int *transa, int *transb, int *m, int *n, int *k, const void *alpha,
   const double *a, int *lda, const double *b, int *ldb, const void *beta,
   void *c, int *ldc) {
  BLAS_zgemm_d_d(blas_colmajor, (enum blas_trans_type) *transa,
		 (enum blas_trans_type) *transb, *m, *n, *k, alpha, a, *lda,
		 b, *ldb, beta, c, *ldc);
}
