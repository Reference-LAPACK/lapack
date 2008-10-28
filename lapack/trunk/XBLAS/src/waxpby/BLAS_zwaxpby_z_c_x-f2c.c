
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby_z_c_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zwaxpby_z_c_x, BLAS_ZWAXPBY_Z_C_X)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw, int *prec) {
  BLAS_zwaxpby_z_c_x(*n, alpha, x, *incx, beta, y, *incy, w, *incw,
		     (enum blas_prec_type) *prec);
}
