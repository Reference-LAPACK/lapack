
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsbmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, double alpha, const float *a, int lda,
		    const float *x, int incx, double beta,
		    double *y, int incy);


extern void FC_FUNC_(blas_dsbmv_s_s, BLAS_DSBMV_S_S)
 
  (int *uplo, int *n, int *k, double *alpha, const float *a, int *lda,
   const float *x, int *incx, double *beta, double *y, int *incy) {
  BLAS_dsbmv_s_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, *alpha,
		 a, *lda, x, *incx, *beta, y, *incy);
}
