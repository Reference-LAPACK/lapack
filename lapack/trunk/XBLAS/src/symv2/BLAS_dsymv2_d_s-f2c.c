
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv2_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const double *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     double beta, double *y, int incy);


extern void FC_FUNC_(blas_dsymv2_d_s, BLAS_DSYMV2_D_S)
 
  (int *uplo, int *n, double *alpha, const double *a, int *lda,
   const float *x_head, const float *x_tail, int *incx, double *beta,
   double *y, int *incy) {
  BLAS_dsymv2_d_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		  *lda, x_head, x_tail, *incx, *beta, y, *incy);
}
