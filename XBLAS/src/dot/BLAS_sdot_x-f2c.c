
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sdot_x(enum blas_conj_type conj, int n, float alpha,
		 const float *x, int incx, float beta,
		 const float *y, int incy,
		 float *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_sdot_x, BLAS_SDOT_X)
 
  (int *conj, int *n, float *alpha, const float *x, int *incx, float *beta,
   const float *y, int *incy, float *r, int *prec) {
  BLAS_sdot_x((enum blas_conj_type) *conj, *n, *alpha, x, *incx, *beta, y,
	      *incy, r, (enum blas_prec_type) *prec);
}
