
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zge_sum_mv_z_c_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const void *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zge_sum_mv_z_c_x, BLAS_ZGE_SUM_MV_Z_C_X)
 
  (int *m, int *n, const void *alpha, const void *a, int *lda, const void *x,
   int *incx, const void *beta, const void *b, int *ldb, void *y, int *incy,
   int *prec) {
  BLAS_zge_sum_mv_z_c_x(blas_colmajor, *m, *n, alpha, a, *lda, x, *incx, beta,
			b, *ldb, y, *incy, (enum blas_prec_type) *prec);
}
