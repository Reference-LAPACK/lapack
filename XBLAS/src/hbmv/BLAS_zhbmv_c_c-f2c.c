
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zhbmv_c_c(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zhbmv_c_c, BLAS_ZHBMV_C_C)
 
  (int *uplo, int *n, int *k, const void *alpha, const void *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zhbmv_c_c(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
