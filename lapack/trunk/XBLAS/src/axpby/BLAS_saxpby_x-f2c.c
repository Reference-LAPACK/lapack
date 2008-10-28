
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_saxpby_x(int n, float alpha, const float *x, int incx,
		   float beta, float *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_saxpby_x, BLAS_SAXPBY_X)
 
  (int *n, float *alpha, const float *x, int *incx, float *beta, float *y,
   int *incy, int *prec) {
  BLAS_saxpby_x(*n, *alpha, x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
