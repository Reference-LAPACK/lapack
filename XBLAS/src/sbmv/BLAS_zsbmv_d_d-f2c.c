
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsbmv_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const double *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zsbmv_d_d, BLAS_ZSBMV_D_D)
 
  (int *uplo, int *n, int *k, const void *alpha, const double *a, int *lda,
   const double *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zsbmv_d_d(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
