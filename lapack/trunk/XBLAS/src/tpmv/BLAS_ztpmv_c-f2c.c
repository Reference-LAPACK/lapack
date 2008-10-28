
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztpmv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp,
		  void *x, int incx);


extern void FC_FUNC_(blas_ztpmv_c, BLAS_ZTPMV_C)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const void *tp, void *x, int *incx) {
  BLAS_ztpmv_c(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       alpha, tp, x, *incx);
}
