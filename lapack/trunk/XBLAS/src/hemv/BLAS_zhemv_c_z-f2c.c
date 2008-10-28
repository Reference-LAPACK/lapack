
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zhemv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zhemv_c_z, BLAS_ZHEMV_C_Z)
 
  (int *uplo, int *n, const void *alpha, const void *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zhemv_c_z(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
