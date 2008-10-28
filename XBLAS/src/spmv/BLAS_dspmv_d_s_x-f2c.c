
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dspmv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const double *ap,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dspmv_d_s_x, BLAS_DSPMV_D_S_X)
 
  (int *uplo, int *n, double *alpha, const double *ap, const float *x,
   int *incx, double *beta, double *y, int *incy, int *prec) {
  BLAS_dspmv_d_s_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap,
		   x, *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
