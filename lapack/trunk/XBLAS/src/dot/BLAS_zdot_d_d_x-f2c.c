
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zdot_d_d_x(enum blas_conj_type conj, int n, const void *alpha,
		     const double *x, int incx, const void *beta,
		     const double *y, int incy,
		     void *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zdot_d_d_x, BLAS_ZDOT_D_D_X)
 
  (int *conj, int *n, const void *alpha, const double *x, int *incx,
   const void *beta, const double *y, int *incy, void *r, int *prec) {
  BLAS_zdot_d_d_x((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		  *incy, r, (enum blas_prec_type) *prec);
}
