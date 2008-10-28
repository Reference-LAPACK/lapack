
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zgemv_d_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const double *a, int lda,
		    const double *x, int incx, const void *beta, void *y,
		    int incy);


extern void FC_FUNC_(blas_zgemv_d_d, BLAS_ZGEMV_D_D)
 
  (int *trans, int *m, int *n, const void *alpha, const double *a, int *lda,
   const double *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zgemv_d_d(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, alpha,
		 a, *lda, x, *incx, beta, y, *incy);
}
