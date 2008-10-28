
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *ap,
		    const float *x, int incx, double beta,
		    double *y, int incy);


extern void FC_FUNC_(blas_dspmv_s_s, BLAS_DSPMV_S_S)
 
  (int *uplo, int *n, double *alpha, const float *ap, const float *x,
   int *incx, double *beta, double *y, int *incy) {
  BLAS_dspmv_s_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap,
		 x, *incx, *beta, y, *incy);
}
