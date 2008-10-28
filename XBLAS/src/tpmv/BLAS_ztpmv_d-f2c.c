
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztpmv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const double *tp,
		  void *x, int incx);


extern void FC_FUNC_(blas_ztpmv_d, BLAS_ZTPMV_D)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const double *tp, void *x, int *incx) {
  BLAS_ztpmv_d(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       alpha, tp, x, *incx);
}
