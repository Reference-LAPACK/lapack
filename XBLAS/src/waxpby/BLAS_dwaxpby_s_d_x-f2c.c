
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dwaxpby_s_d_x(int n, double alpha, const float *x, int incx,
			double beta, const double *y, int incy, double *w,
			int incw, enum blas_prec_type prec);


extern void FC_FUNC_(blas_dwaxpby_s_d_x, BLAS_DWAXPBY_S_D_X)
 
  (int *n, double *alpha, const float *x, int *incx, double *beta,
   const double *y, int *incy, double *w, int *incw, int *prec) {
  BLAS_dwaxpby_s_d_x(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw,
		     (enum blas_prec_type) *prec);
}
