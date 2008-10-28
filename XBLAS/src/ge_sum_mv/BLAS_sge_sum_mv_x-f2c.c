
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sge_sum_mv_x(enum blas_order_type order, int m, int n,
		       float alpha, const float *a, int lda,
		       const float *x, int incx,
		       float beta, const float *b, int ldb,
		       float *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_sge_sum_mv_x, BLAS_SGE_SUM_MV_X)
 
  (int *m, int *n, float *alpha, const float *a, int *lda, const float *x,
   int *incx, float *beta, const float *b, int *ldb, float *y, int *incy,
   int *prec) {
  BLAS_sge_sum_mv_x(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta,
		    b, *ldb, y, *incy, (enum blas_prec_type) *prec);
}
