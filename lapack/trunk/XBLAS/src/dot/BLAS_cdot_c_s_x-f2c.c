
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cdot_c_s_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const float *y, int incy,
		     void *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cdot_c_s_x, BLAS_CDOT_C_S_X)
 
  (int *conj, int *n, const void *alpha, const void *x, int *incx,
   const void *beta, const float *y, int *incy, void *r, int *prec) {
  BLAS_cdot_c_s_x((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		  *incy, r, (enum blas_prec_type) *prec);
}
