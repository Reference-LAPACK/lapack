
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ddot_s_s(enum blas_conj_type conj, int n, double alpha,
		   const float *x, int incx, double beta,
		   const float *y, int incy, double *r);


extern void FC_FUNC_(blas_ddot_s_s, BLAS_DDOT_S_S)
 
  (int *conj, int *n, double *alpha, const float *x, int *incx, double *beta,
   const float *y, int *incy, double *r) {
  BLAS_ddot_s_s((enum blas_conj_type) *conj, *n, *alpha, x, *incx, *beta, y,
		*incy, r);
}
