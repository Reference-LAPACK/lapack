#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zaxpby_c(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y, int incy)

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
 * alpha     (input) const void*
 *
 * x         (input) const void*
 *           Array of length n.
 *
 * incx      (input) int
 *           The stride used to access components x[i].
 * 
 * beta      (input) const void*
 *
 * y         (input) void*
 *           Array of length n.
 * 
 * incy      (input) int
 *           The stride used to access components y[i].
 *
 */
{
  static const char routine_name[] = "BLAS_zaxpby_c";

  int i, ix = 0, iy = 0;
  const float *x_i = (float *) x;
  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float x_ii[2];
  double y_ii[2];
  double tmpx[2];
  double tmpy[2];


  /* Test the input parameters. */
  if (incx == 0)
    BLAS_error(routine_name, -4, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -7, incy, NULL);

  /* Immediate return */
  if (n <= 0
      || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)))
    return;



  incx *= 2;
  incy *= 2;
  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = y_i[iy];
    y_ii[1] = y_i[iy + 1];
    {
      tmpx[0] = (double) alpha_i[0] * x_ii[0] - (double) alpha_i[1] * x_ii[1];
      tmpx[1] = (double) alpha_i[0] * x_ii[1] + (double) alpha_i[1] * x_ii[0];
    }				/* tmpx  = alpha * x[ix] */
    {
      tmpy[0] = (double) beta_i[0] * y_ii[0] - (double) beta_i[1] * y_ii[1];
      tmpy[1] = (double) beta_i[0] * y_ii[1] + (double) beta_i[1] * y_ii[0];
    }				/* tmpy = beta * y[iy] */
    tmpy[0] = tmpy[0] + tmpx[0];
    tmpy[1] = tmpy[1] + tmpx[1];
    y_i[iy] = tmpy[0];
    y_i[iy + 1] = tmpy[1];
    ix += incx;
    iy += incy;
  }				/* endfor */



}				/* end BLAS_zaxpby_c */
