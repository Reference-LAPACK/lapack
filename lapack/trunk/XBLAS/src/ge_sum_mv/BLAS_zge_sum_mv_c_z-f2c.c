
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zge_sum_mv_c_z(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const void *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);


extern void FC_FUNC_(blas_zge_sum_mv_c_z, BLAS_ZGE_SUM_MV_C_Z)
 
  (int *m, int *n, const void *alpha, const void *a, int *lda, const void *x,
   int *incx, const void *beta, const void *b, int *ldb, void *y, int *incy) {
  BLAS_zge_sum_mv_c_z(blas_colmajor, *m, *n, alpha, a, *lda, x, *incx, beta,
		      b, *ldb, y, *incy);
}
