#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_zwaxpby_d_d(int n, const void *alpha, const double *x, int incx,
		      const void *beta, const double *y, int incy, void *w,
		      int incw)

/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *     w <- alpha * x + beta * y
 * 
 * Arguments
 * =========
 *
 * n     (input) int
 *       The length of vectors x, y, and w.
 * 
 * alpha (input) const void*
 *
 * x     (input) const double*
 *       Array of length n.
 * 
 * incx  (input) int
 *       The stride used to access components x[i].
 *
 * beta  (input) const void*
 *
 * y     (input) double*
 *       Array of length n.
 *
 * incy  (input) int
 *       The stride used to access components y[i].
 *
 * w     (output) void*
 *       Array of length n.
 *
 * incw  (input) int
 *       The stride used to write components w[i].
 *
 */
{
  char *routine_name = "BLAS_zwaxpby_d_d";

  int i, ix = 0, iy = 0, iw = 0;
  double *w_i = (double *) w;
  const double *x_i = x;
  const double *y_i = y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double x_ii;
  double y_ii;
  double tmpx[2];
  double tmpy[2];



  /* Test the input parameters. */
  if (incx == 0)
    BLAS_error(routine_name, -4, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -7, incy, NULL);
  else if (incw == 0)
    BLAS_error(routine_name, -9, incw, NULL);


  /* Immediate return */
  if (n <= 0) {
    return;
  }





  incw *= 2;
  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;
  if (incw < 0)
    iw = (-n + 1) * incw;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];
    {
      tmpx[0] = alpha_i[0] * x_ii;
      tmpx[1] = alpha_i[1] * x_ii;
    }				/* tmpx  = alpha * x[ix] */
    {
      tmpy[0] = beta_i[0] * y_ii;
      tmpy[1] = beta_i[1] * y_ii;
    }				/* tmpy = beta * y[iy] */
    tmpy[0] = tmpy[0] + tmpx[0];
    tmpy[1] = tmpy[1] + tmpx[1];
    w_i[iw] = tmpy[0];
    w_i[iw + 1] = tmpy[1];
    ix += incx;
    iy += incy;
    iw += incw;
  }				/* endfor */



}
