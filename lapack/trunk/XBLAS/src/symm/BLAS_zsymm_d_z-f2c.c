
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymm_d_z(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const double *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);


extern void FC_FUNC_(blas_zsymm_d_z, BLAS_ZSYMM_D_Z)
 
  (int *side, int *uplo, int *m, int *n, const void *alpha, const double *a,
   int *lda, const void *b, int *ldb, const void *beta, void *c, int *ldc) {
  BLAS_zsymm_d_z(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, alpha, a, *lda, b, *ldb,
		 beta, c, *ldc);
}
