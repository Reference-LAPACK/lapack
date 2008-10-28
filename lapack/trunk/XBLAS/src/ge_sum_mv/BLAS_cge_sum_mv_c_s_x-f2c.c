
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cge_sum_mv_c_s_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const float *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cge_sum_mv_c_s_x, BLAS_CGE_SUM_MV_C_S_X)
 
  (int *m, int *n, const void *alpha, const void *a, int *lda, const float *x,
   int *incx, const void *beta, const void *b, int *ldb, void *y, int *incy,
   int *prec) {
  BLAS_cge_sum_mv_c_s_x(blas_colmajor, *m, *n, alpha, a, *lda, x, *incx, beta,
			b, *ldb, y, *incy, (enum blas_prec_type) *prec);
}
