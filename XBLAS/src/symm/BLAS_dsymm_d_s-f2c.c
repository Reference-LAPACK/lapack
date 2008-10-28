
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymm_d_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const double *a, int lda,
		    const float *b, int ldb, double beta, double *c, int ldc);


extern void FC_FUNC_(blas_dsymm_d_s, BLAS_DSYMM_D_S)
 
  (int *side, int *uplo, int *m, int *n, double *alpha, const double *a,
   int *lda, const float *b, int *ldb, double *beta, double *c, int *ldc) {
  BLAS_dsymm_d_s(blas_colmajor, (enum blas_side_type) *side,
		 (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b,
		 *ldb, *beta, c, *ldc);
}
