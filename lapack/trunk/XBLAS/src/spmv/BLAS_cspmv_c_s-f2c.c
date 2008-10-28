
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cspmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_cspmv_c_s, BLAS_CSPMV_C_S)
 
  (int *uplo, int *n, const void *alpha, const void *ap, const float *x,
   int *incx, const void *beta, void *y, int *incy) {
  BLAS_cspmv_c_s(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, ap, x,
		 *incx, beta, y, *incy);
}
