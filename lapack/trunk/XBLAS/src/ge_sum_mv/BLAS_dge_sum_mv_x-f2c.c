
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dge_sum_mv_x(enum blas_order_type order, int m, int n,
		       double alpha, const double *a, int lda,
		       const double *x, int incx,
		       double beta, const double *b, int ldb,
		       double *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dge_sum_mv_x, BLAS_DGE_SUM_MV_X)
 
  (int *m, int *n, double *alpha, const double *a, int *lda, const double *x,
   int *incx, double *beta, const double *b, int *ldb, double *y, int *incy,
   int *prec) {
  BLAS_dge_sum_mv_x(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta,
		    b, *ldb, y, *incy, (enum blas_prec_type) *prec);
}
