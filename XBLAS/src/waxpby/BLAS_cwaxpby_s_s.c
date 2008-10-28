#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_cwaxpby_s_s(int n, const void *alpha, const float *x, int incx,
		      const void *beta, const float *y, int incy, void *w,
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
 * x     (input) const float*
 *       Array of length n.
 * 
 * incx  (input) int
 *       The stride used to access components x[i].
 *
 * beta  (input) const void*
 *
 * y     (input) float*
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
  char *routine_name = "BLAS_cwaxpby_s_s";

  int i, ix = 0, iy = 0, iw = 0;
  float *w_i = (float *) w;
  const float *x_i = x;
  const float *y_i = y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float x_ii;
  float y_ii;
  float tmpx[2];
  float tmpy[2];



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
