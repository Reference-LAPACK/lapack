
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dtbsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, double alpha, const float *t, int ldt,
		  double *x, int incx);


extern void FC_FUNC_(blas_dtbsv_s, BLAS_DTBSV_S)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, double *alpha,
   const float *t, int *ldt, double *x, int *incx) {
  BLAS_dtbsv_s(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *k, *alpha, t, *ldt, x, *incx);
}
