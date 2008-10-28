
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ctpmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const float *tp,
		  void *x, int incx);


extern void FC_FUNC_(blas_ctpmv_s, BLAS_CTPMV_S)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const float *tp, void *x, int *incx) {
  BLAS_ctpmv_s(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       alpha, tp, x, *incx);
}
