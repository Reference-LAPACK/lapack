
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztbsv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const double *t, int ldt,
		  void *x, int incx);


extern void FC_FUNC_(blas_ztbsv_d, BLAS_ZTBSV_D)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, const void *alpha,
   const double *t, int *ldt, void *x, int *incx) {
  BLAS_ztbsv_d(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *k, alpha, t, *ldt, x, *incx);
}
