
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztrsv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *T, int ldt,
		  void *x, int incx);


extern void FC_FUNC_(blas_ztrsv_c, BLAS_ZTRSV_C)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha, const void *T,
   int *ldt, void *x, int *incx) {
  BLAS_ztrsv_c(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       alpha, T, *ldt, x, *incx);
}
