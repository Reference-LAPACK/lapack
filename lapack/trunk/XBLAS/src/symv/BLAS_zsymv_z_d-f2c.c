
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zsymv_z_d, BLAS_ZSYMV_Z_D)
 
  (int *uplo, int *n, const void *alpha, const void *a, int *lda,
   const double *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zsymv_z_d(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
