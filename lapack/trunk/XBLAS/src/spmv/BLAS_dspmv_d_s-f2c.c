
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dspmv_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const double *ap,
		    const float *x, int incx, double beta,
		    double *y, int incy);


extern void FC_FUNC_(blas_dspmv_d_s, BLAS_DSPMV_D_S)
 
  (int *uplo, int *n, double *alpha, const double *ap, const float *x,
   int *incx, double *beta, double *y, int *incy) {
  BLAS_dspmv_d_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap,
		 x, *incx, *beta, y, *incy);
}
