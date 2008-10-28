
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dge_sum_mv_s_d(enum blas_order_type order, int m, int n,
			 double alpha, const float *a, int lda,
			 const double *x, int incx,
			 double beta, const float *b, int ldb,
			 double *y, int incy);


extern void FC_FUNC_(blas_dge_sum_mv_s_d, BLAS_DGE_SUM_MV_S_D)
 
  (int *m, int *n, double *alpha, const float *a, int *lda, const double *x,
   int *incx, double *beta, const float *b, int *ldb, double *y, int *incy) {
  BLAS_dge_sum_mv_s_d(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta,
		      b, *ldb, y, *incy);
}
