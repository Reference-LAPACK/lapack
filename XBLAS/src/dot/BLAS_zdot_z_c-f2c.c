
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zdot_z_c(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const void *y, int incy, void *r);


extern void FC_FUNC_(blas_zdot_z_c, BLAS_ZDOT_Z_C)
 
  (int *conj, int *n, const void *alpha, const void *x, int *incx,
   const void *beta, const void *y, int *incy, void *r) {
  BLAS_zdot_z_c((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		*incy, r);
}
