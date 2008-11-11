
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csymv2_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const float *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_csymv2_s_c, BLAS_CSYMV2_S_C)
 
  (int *uplo, int *n, const void *alpha, const float *a, int *lda,
   const void *x_head, const void *x_tail, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_csymv2_s_c(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		  *lda, x_head, x_tail, *incx, beta, y, *incy);
}
