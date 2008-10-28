
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zge_sum_mv_d_d(enum blas_order_type order, int m, int n,
			 const void *alpha, const double *a, int lda,
			 const double *x, int incx,
			 const void *beta, const double *b, int ldb,
			 void *y, int incy);


extern void FC_FUNC_(blas_zge_sum_mv_d_d, BLAS_ZGE_SUM_MV_D_D)
 
  (int *m, int *n, const void *alpha, const double *a, int *lda,
   const double *x, int *incx, const void *beta, const double *b, int *ldb,
   void *y, int *incy) {
  BLAS_zge_sum_mv_d_d(blas_colmajor, *m, *n, alpha, a, *lda, x, *incx, beta,
		      b, *ldb, y, *incy);
}
