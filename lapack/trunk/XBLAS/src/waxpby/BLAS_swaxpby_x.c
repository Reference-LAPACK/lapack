#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_swaxpby_x(int n, float alpha, const float *x, int incx,
		    float beta, const float *y, int incy, float *w,
		    int incw, enum blas_prec_type prec)

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
 * alpha (input) float
 *
 * x     (input) const float*
 *       Array of length n.
 * 
 * incx  (input) int
 *       The stride used to access components x[i].
 *
 * beta  (input) float
 *
 * y     (input) float*
 *       Array of length n.
 *
 * incy  (input) int
 *       The stride used to access components y[i].
 *
 * w     (output) float*
 *       Array of length n.
 *
 * incw  (input) int
 *       The stride used to write components w[i].
 *
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  char *routine_name = "BLAS_swaxpby_x";
  switch (prec) {
  case blas_prec_single:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = w;
      const float *x_i = x;
      const float *y_i = y;
      float alpha_i = alpha;
      float beta_i = beta;
      float x_ii;
      float y_ii;
      float tmpx;
      float tmpy;



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






      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];
	tmpx = alpha_i * x_ii;	/* tmpx  = alpha * x[ix] */
	tmpy = beta_i * y_ii;	/* tmpy = beta * y[iy] */
	tmpy = tmpy + tmpx;
	w_i[iw] = tmpy;
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */



      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = w;
      const float *x_i = x;
      const float *y_i = y;
      float alpha_i = alpha;
      float beta_i = beta;
      float x_ii;
      float y_ii;
      double tmpx;
      double tmpy;



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






      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];
	tmpx = (double) alpha_i *x_ii;	/* tmpx  = alpha * x[ix] */
	tmpy = (double) beta_i *y_ii;	/* tmpy = beta * y[iy] */
	tmpy = tmpy + tmpx;
	w_i[iw] = tmpy;
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */



      break;
    }

  case blas_prec_extra:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = w;
      const float *x_i = x;
      const float *y_i = y;
      float alpha_i = alpha;
      float beta_i = beta;
      float x_ii;
      float y_ii;
      double head_tmpx, tail_tmpx;
      double head_tmpy, tail_tmpy;

      FPU_FIX_DECL;

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

      FPU_FIX_START;




      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];
	head_tmpx = (double) alpha_i *x_ii;
	tail_tmpx = 0.0;	/* tmpx  = alpha * x[ix] */
	head_tmpy = (double) beta_i *y_ii;
	tail_tmpy = 0.0;	/* tmpy = beta * y[iy] */
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_tmpy + head_tmpx;
	  bv = s1 - head_tmpy;
	  s2 = ((head_tmpx - bv) + (head_tmpy - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_tmpy + tail_tmpx;
	  bv = t1 - tail_tmpy;
	  t2 = ((tail_tmpx - bv) + (tail_tmpy - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_tmpy = t1 + t2;
	  tail_tmpy = t2 - (head_tmpy - t1);
	}
	w_i[iw] = head_tmpy;
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */

      FPU_FIX_STOP;

      break;
    }
  }
}
