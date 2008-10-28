
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsbmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zsbmv_c_z, BLAS_ZSBMV_C_Z)
 
  (int *uplo, int *n, int *k, const void *alpha, const void *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zsbmv_c_z(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
