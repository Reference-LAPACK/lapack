
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zhpmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);


extern void FC_FUNC_(blas_zhpmv_c_z, BLAS_ZHPMV_C_Z)
 
  (int *uplo, int *n, const void *alpha, const void *ap, const void *x,
   int *incx, const void *beta, void *y, int *incy) {
  BLAS_zhpmv_c_z(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, ap, x,
		 *incx, beta, y, *incy);
}
