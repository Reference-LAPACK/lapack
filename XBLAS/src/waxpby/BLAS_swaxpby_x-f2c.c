
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_swaxpby_x(int n, float alpha, const float *x, int incx,
		    float beta, const float *y, int incy, float *w,
		    int incw, enum blas_prec_type prec);


extern void FC_FUNC_(blas_swaxpby_x, BLAS_SWAXPBY_X)
 
  (int *n, float *alpha, const float *x, int *incx, float *beta,
   const float *y, int *incy, float *w, int *incw, int *prec) {
  BLAS_swaxpby_x(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw,
		 (enum blas_prec_type) *prec);
}
