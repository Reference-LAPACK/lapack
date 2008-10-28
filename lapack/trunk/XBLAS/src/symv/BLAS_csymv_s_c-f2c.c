
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csymv_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_csymv_s_c, BLAS_CSYMV_S_C)
 
  (int *uplo, int *n, const void *alpha, const float *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_csymv_s_c(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
