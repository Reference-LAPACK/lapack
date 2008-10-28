
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymm_s_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const float *a, int lda,
		    const float *b, int ldb, double beta, double *c, int ldc);


extern void FC_FUNC_(blas_dsymm_s_s, BLAS_DSYMM_S_S)
 
  (int *side, int *uplo, int *m, int *n, double *alpha, const float *a,
   int *lda, const float *b, int *ldb, double *beta, double *c, int *ldc) {
  BLAS_dsymm_s_s(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b,
		 *ldb, *beta, c, *ldc);
}
