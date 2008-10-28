
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dtrmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, double alpha, const float *T, int ldt, double *x,
		    int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dtrmv_s_x, BLAS_DTRMV_S_X)
 
  (int *uplo, int *trans, int *diag, int *n, double *alpha, const float *T,
   int *ldt, double *x, int *incx, int *prec) {
  BLAS_dtrmv_s_x(blas_colmajor, (enum blas_uplo_type) *uplo,
		 (enum blas_trans_type) *trans, (enum blas_diag_type) *diag,
		 *n, *alpha, T, *ldt, x, *incx, (enum blas_prec_type) *prec);
}
