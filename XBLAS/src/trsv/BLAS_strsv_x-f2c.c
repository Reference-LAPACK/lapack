
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_strsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, float alpha, const float *T, int ldt,
		  float *x, int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_strsv_x, BLAS_STRSV_X)
 
  (int *uplo, int *trans, int *diag, int *n, float *alpha, const float *T,
   int *ldt, float *x, int *incx, int *prec) {
  BLAS_strsv_x(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *alpha, T, *ldt, x, *incx, (enum blas_prec_type) *prec);
}
