
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta,
		    double *y, int incy);


extern void FC_FUNC_(blas_dsymv_s_d, BLAS_DSYMV_S_D)
 
  (int *uplo, int *n, double *alpha, const float *a, int *lda,
   const double *x, int *incx, double *beta, double *y, int *incy) {
  BLAS_dsymv_s_d(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		 *lda, x, *incx, *beta, y, *incy);
}
