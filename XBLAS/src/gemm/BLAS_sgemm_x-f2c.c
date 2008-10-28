
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  float alpha, const float *a, int lda, const float *b,
		  int ldb, float beta, float *c, int ldc,
		  enum blas_prec_type prec);


extern void FC_FUNC_(blas_sgemm_x, BLAS_SGEMM_X)
 
  (int *transa, int *transb, int *m, int *n, int *k, float *alpha,
   const float *a, int *lda, const float *b, int *ldb, float *beta, float *c,
   int *ldc, int *prec) {
  BLAS_sgemm_x(blas_colmajor, (enum blas_trans_type) *transa,
	       (enum blas_trans_type) *transb, *m, *n, *k, *alpha, a, *lda, b,
	       *ldb, *beta, c, *ldc, (enum blas_prec_type) *prec);
}
