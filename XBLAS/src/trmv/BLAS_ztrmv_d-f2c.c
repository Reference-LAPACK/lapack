
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztrmv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const double *T, int ldt,
		  void *x, int incx);


extern void FC_FUNC_(blas_ztrmv_d, BLAS_ZTRMV_D)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const double *T, int *ldt, void *x, int *incx) {
  BLAS_ztrmv_d(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       alpha, T, *ldt, x, *incx);
}
