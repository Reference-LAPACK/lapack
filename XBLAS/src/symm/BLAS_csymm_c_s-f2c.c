
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csymm_c_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const float *b, int ldb, const void *beta,
		    void *c, int ldc);


extern void FC_FUNC_(blas_csymm_c_s, BLAS_CSYMM_C_S)
 
  (int *side, int *uplo, int *m, int *n, const void *alpha, const void *a,
   int *lda, const float *b, int *ldb, const void *beta, void *c, int *ldc) {
  BLAS_csymm_c_s(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, alpha, a, *lda, b, *ldb,
		 beta, c, *ldc);
}
