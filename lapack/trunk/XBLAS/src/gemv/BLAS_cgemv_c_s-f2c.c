
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemv_c_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta, void *y,
		    int incy);


extern void FC_FUNC_(blas_cgemv_c_s, BLAS_CGEMV_C_S)
 
  (int *trans, int *m, int *n, const void *alpha, const void *a, int *lda,
   const float *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_cgemv_c_s(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, alpha,
		 a, *lda, x, *incx, beta, y, *incy);
}
