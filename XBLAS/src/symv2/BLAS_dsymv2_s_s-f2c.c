
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv2_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const float *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     double beta, double *y, int incy);


extern void FC_FUNC_(blas_dsymv2_s_s, BLAS_DSYMV2_S_S)
 
  (int *uplo, int *n, double *alpha, const float *a, int *lda,
   const float *x_head, const float *x_tail, int *incx, double *beta,
   double *y, int *incy) {
  BLAS_dsymv2_s_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		  *lda, x_head, x_tail, *incx, *beta, y, *incy);
}
