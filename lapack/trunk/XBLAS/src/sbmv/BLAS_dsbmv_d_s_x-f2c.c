
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsbmv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, double alpha, const double *a, int lda,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dsbmv_d_s_x, BLAS_DSBMV_D_S_X)
 
  (int *uplo, int *n, int *k, double *alpha, const double *a, int *lda,
   const float *x, int *incx, double *beta, double *y, int *incy, int *prec) {
  BLAS_dsbmv_d_s_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, *alpha,
		   a, *lda, x, *incx, *beta, y, *incy,
		   (enum blas_prec_type) *prec);
}
