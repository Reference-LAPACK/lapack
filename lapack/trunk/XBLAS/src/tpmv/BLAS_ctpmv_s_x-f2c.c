
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ctpmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const float *tp,
		    void *x, int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ctpmv_s_x, BLAS_CTPMV_S_X)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const float *tp, void *x, int *incx, int *prec) {
  BLAS_ctpmv_s_x(blas_colmajor, (enum blas_uplo_type) *uplo,
		 (enum blas_trans_type) *trans, (enum blas_diag_type) *diag,
		 *n, alpha, tp, x, *incx, (enum blas_prec_type) *prec);
}
