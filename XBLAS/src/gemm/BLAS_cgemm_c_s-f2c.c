
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemm_c_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda, const float *b,
		    int ldb, const void *beta, void *c, int ldc);


extern void FC_FUNC_(blas_cgemm_c_s, BLAS_CGEMM_C_S)
 
  (int *transa, int *transb, int *m, int *n, int *k, const void *alpha,
   const void *a, int *lda, const float *b, int *ldb, const void *beta,
   void *c, int *ldc) {
  BLAS_cgemm_c_s(blas_colmajor, (enum blas_trans_type) *transa,
		 (enum blas_trans_type) *transb, *m, *n, *k, alpha, a, *lda,
		 b, *ldb, beta, c, *ldc);
}
