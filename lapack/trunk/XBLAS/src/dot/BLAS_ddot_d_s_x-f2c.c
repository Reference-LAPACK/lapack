
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ddot_d_s_x(enum blas_conj_type conj, int n, double alpha,
		     const double *x, int incx, double beta,
		     const float *y, int incy,
		     double *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ddot_d_s_x, BLAS_DDOT_D_S_X)
 
  (int *conj, int *n, double *alpha, const double *x, int *incx, double *beta,
   const float *y, int *incy, double *r, int *prec) {
  BLAS_ddot_d_s_x((enum blas_conj_type) *conj, *n, *alpha, x, *incx, *beta, y,
		  *incy, r, (enum blas_prec_type) *prec);
}
