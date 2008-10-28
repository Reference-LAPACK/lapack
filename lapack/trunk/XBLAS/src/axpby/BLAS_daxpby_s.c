#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_daxpby_s(int n, double alpha, const float *x, int incx,
		   double beta, double *y, int incy)

/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *      y <- alpha * x + beta * y.
 *
 * Arguments
 * =========
 * 
 * n         (input) int
 *           The length of vectors x and y.
 * 
 * alpha     (input) double
 *
 * x         (input) const float*
 *           Array of length n.
 *
 * incx      (input) int
 *           The stride used to access components x[i].
 * 
 * beta      (input) double
 *
 * y         (input) double*
 *           Array of length n.
 * 
 * incy      (input) int
 *           The stride used to access components y[i].
 *
 */
{
  static const char routine_name[] = "BLAS_daxpby_s";

  int i, ix = 0, iy = 0;
  const float *x_i = x;
  double *y_i = y;
  double alpha_i = alpha;
  double beta_i = beta;
  float x_ii;
  double y_ii;
  double tmpx;
  double tmpy;


  /* Test the input parameters. */
  if (incx == 0)
    BLAS_error(routine_name, -4, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -7, incy, NULL);

  /* Immediate return */
  if (n <= 0 || (alpha_i == 0.0 && beta_i == 1.0))
    return;





  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];
    tmpx = alpha_i * x_ii;	/* tmpx  = alpha * x[ix] */
    tmpy = beta_i * y_ii;	/* tmpy = beta * y[iy] */
    tmpy = tmpy + tmpx;
    y_i[iy] = tmpy;
    ix += incx;
    iy += incy;
  }				/* endfor */



}				/* end BLAS_daxpby_s */
