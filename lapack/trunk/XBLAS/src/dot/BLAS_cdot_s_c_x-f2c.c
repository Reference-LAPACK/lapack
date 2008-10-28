
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cdot_s_c_x(enum blas_conj_type conj, int n, const void *alpha,
		     const float *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cdot_s_c_x, BLAS_CDOT_S_C_X)
 
  (int *conj, int *n, const void *alpha, const float *x, int *incx,
   const void *beta, const void *y, int *incy, void *r, int *prec) {
  BLAS_cdot_s_c_x((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		  *incy, r, (enum blas_prec_type) *prec);
}
