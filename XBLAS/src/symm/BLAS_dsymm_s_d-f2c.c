
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymm_s_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const float *a, int lda,
		    const double *b, int ldb, double beta,
		    double *c, int ldc);


extern void FC_FUNC_(blas_dsymm_s_d, BLAS_DSYMM_S_D)
 
  (int *side, int *uplo, int *m, int *n, double *alpha, const float *a,
   int *lda, const double *b, int *ldb, double *beta, double *c, int *ldc) {
  BLAS_dsymm_s_d(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b,
		 *ldb, *beta, c, *ldc);
}
