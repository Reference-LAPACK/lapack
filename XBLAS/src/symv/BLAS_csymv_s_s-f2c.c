
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csymv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_csymv_s_s, BLAS_CSYMV_S_S)
 
  (int *uplo, int *n, const void *alpha, const float *a, int *lda,
   const float *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_csymv_s_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
