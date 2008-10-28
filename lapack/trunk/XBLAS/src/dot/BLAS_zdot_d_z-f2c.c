
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zdot_d_z(enum blas_conj_type conj, int n, const void *alpha,
		   const double *x, int incx, const void *beta,
		   const void *y, int incy, void *r);


extern void FC_FUNC_(blas_zdot_d_z, BLAS_ZDOT_D_Z)
 
  (int *conj, int *n, const void *alpha, const double *x, int *incx,
   const void *beta, const void *y, int *incy, void *r) {
  BLAS_zdot_d_z((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		*incy, r);
}
