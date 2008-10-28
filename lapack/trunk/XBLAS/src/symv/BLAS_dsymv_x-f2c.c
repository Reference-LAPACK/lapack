
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, double alpha, const double *a, int lda,
		  const double *x, int incx, double beta,
		  double *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dsymv_x, BLAS_DSYMV_X)
 
  (int *uplo, int *n, double *alpha, const double *a, int *lda,
   const double *x, int *incx, double *beta, double *y, int *incy,
   int *prec) {
  BLAS_dsymv_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
	       *lda, x, *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
