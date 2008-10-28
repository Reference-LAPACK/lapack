
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zhemm_z_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const double *b, int ldb, const void *beta,
		    void *c, int ldc);


extern void FC_FUNC_(blas_zhemm_z_d, BLAS_ZHEMM_Z_D)
 
  (int *side, int *uplo, int *m, int *n, const void *alpha, const void *a,
   int *lda, const double *b, int *ldb, const void *beta, void *c, int *ldc) {
  BLAS_zhemm_z_d(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, alpha, a, *lda, b, *ldb,
		 beta, c, *ldc);
}
