
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_chemv2_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     const void *beta, const float *y, int incy);


extern void FC_FUNC_(blas_chemv2_c_s, BLAS_CHEMV2_C_S)
 
  (int *uplo, int *n, const void *alpha, const void *a, int *lda,
   const float *x_head, const float *x_tail, int *incx, const void *beta,
   const float *y, int *incy) {
  BLAS_chemv2_c_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		  *lda, x_head, x_tail, *incx, beta, y, *incy);
}
