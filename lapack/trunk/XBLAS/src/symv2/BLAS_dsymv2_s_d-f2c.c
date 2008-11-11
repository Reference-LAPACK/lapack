
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv2_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const float *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     double beta, double *y, int incy);


extern void FC_FUNC_(blas_dsymv2_s_d, BLAS_DSYMV2_S_D)
 
  (int *uplo, int *n, double *alpha, const float *a, int *lda,
   const double *x_head, const double *x_tail, int *incx, double *beta,
   double *y, int *incy) {
  BLAS_dsymv2_s_d(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		  *lda, x_head, x_tail, *incx, *beta, y, *incy);
}
