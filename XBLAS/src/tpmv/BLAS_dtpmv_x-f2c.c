
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dtpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const double *tp,
		  double *x, int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dtpmv_x, BLAS_DTPMV_X)
 
  (int *uplo, int *trans, int *diag, int *n, double *alpha, const double *tp,
   double *x, int *incx, int *prec) {
  BLAS_dtpmv_x(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *alpha, tp, x, *incx, (enum blas_prec_type) *prec);
}
