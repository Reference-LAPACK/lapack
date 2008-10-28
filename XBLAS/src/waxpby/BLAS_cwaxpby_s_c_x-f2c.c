
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cwaxpby_s_c_x(int n, const void *alpha, const float *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cwaxpby_s_c_x, BLAS_CWAXPBY_S_C_X)
 
  (int *n, const void *alpha, const float *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw, int *prec) {
  BLAS_cwaxpby_s_c_x(*n, alpha, x, *incx, beta, y, *incy, w, *incw,
		     (enum blas_prec_type) *prec);
}
