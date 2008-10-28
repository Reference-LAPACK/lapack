
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymm_s_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      double alpha, const float *a, int lda,
		      const float *b, int ldb, double beta,
		      double *c, int ldc, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dsymm_s_s_x, BLAS_DSYMM_S_S_X)
 
  (int *side, int *uplo, int *m, int *n, double *alpha, const float *a,
   int *lda, const float *b, int *ldb, double *beta, double *c, int *ldc,
   int *prec) {
  BLAS_dsymm_s_s_x(blas_colmajor, (enum blas_side_type) *side,
		   (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b,
		   *ldb, *beta, c, *ldc, (enum blas_prec_type) *prec);
}
