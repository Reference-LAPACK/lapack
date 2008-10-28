
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dspmv_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const float *ap,
		      const double *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dspmv_s_d_x, BLAS_DSPMV_S_D_X)
 
  (int *uplo, int *n, double *alpha, const float *ap, const double *x,
   int *incx, double *beta, double *y, int *incy, int *prec) {
  BLAS_dspmv_s_d_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap,
		   x, *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
