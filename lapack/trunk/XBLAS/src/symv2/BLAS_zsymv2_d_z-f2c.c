
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymv2_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const double *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_zsymv2_d_z, BLAS_ZSYMV2_D_Z)
 
  (int *uplo, int *n, const void *alpha, const double *a, int *lda,
   const void *x_head, const void *x_tail, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_zsymv2_d_z(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		  *lda, x_head, x_tail, *incx, beta, y, *incy);
}
