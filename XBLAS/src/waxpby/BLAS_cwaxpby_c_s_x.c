#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_cwaxpby_c_s_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const float *y, int incy, void *w,
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
 * alpha (input) const void*
 *
 * x     (input) const void*
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
  char *routine_name = "BLAS_cwaxpby_c_s_x";
  switch (prec) {
  case blas_prec_single:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = (float *) w;
      const float *x_i = (float *) x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii[2];
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



      incx *= 2;

      incw *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii[0] = x_i[ix];
	x_ii[1] = x_i[ix + 1];
	y_ii = y_i[iy];
	{
	  tmpx[0] = alpha_i[0] * x_ii[0] - alpha_i[1] * x_ii[1];
	  tmpx[1] = alpha_i[0] * x_ii[1] + alpha_i[1] * x_ii[0];
	}
	/* tmpx  = alpha * x[ix] */
	{
	  tmpy[0] = beta_i[0] * y_ii;
	  tmpy[1] = beta_i[1] * y_ii;
	}			/* tmpy = beta * y[iy] */
	tmpy[0] = tmpy[0] + tmpx[0];
	tmpy[1] = tmpy[1] + tmpx[1];
	w_i[iw] = tmpy[0];
	w_i[iw + 1] = tmpy[1];
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */



      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = (float *) w;
      const float *x_i = (float *) x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii[2];
      float y_ii;
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



      incx *= 2;

      incw *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii[0] = x_i[ix];
	x_ii[1] = x_i[ix + 1];
	y_ii = y_i[iy];
	{
	  tmpx[0] =
	    (double) alpha_i[0] * x_ii[0] - (double) alpha_i[1] * x_ii[1];
	  tmpx[1] =
	    (double) alpha_i[0] * x_ii[1] + (double) alpha_i[1] * x_ii[0];
	}			/* tmpx  = alpha * x[ix] */
	{
	  tmpy[0] = (double) beta_i[0] * y_ii;
	  tmpy[1] = (double) beta_i[1] * y_ii;
	}			/* tmpy = beta * y[iy] */
	tmpy[0] = tmpy[0] + tmpx[0];
	tmpy[1] = tmpy[1] + tmpx[1];
	w_i[iw] = tmpy[0];
	w_i[iw + 1] = tmpy[1];
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */



      break;
    }

  case blas_prec_extra:{

      int i, ix = 0, iy = 0, iw = 0;
      float *w_i = (float *) w;
      const float *x_i = (float *) x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii[2];
      float y_ii;
      double head_tmpx[2], tail_tmpx[2];
      double head_tmpy[2], tail_tmpy[2];

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

      incx *= 2;

      incw *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;
      if (incw < 0)
	iw = (-n + 1) * incw;

      for (i = 0; i < n; ++i) {
	x_ii[0] = x_i[ix];
	x_ii[1] = x_i[ix + 1];
	y_ii = y_i[iy];
	{
	  double head_e1, tail_e1;
	  double d1;
	  double d2;
	  /* Real part */
	  d1 = (double) alpha_i[0] * x_ii[0];
	  d2 = (double) -alpha_i[1] * x_ii[1];
	  {
	    /* Compute double-double = double + double. */
	    double e, t1, t2;

	    /* Knuth trick. */
	    t1 = d1 + d2;
	    e = t1 - d1;
	    t2 = ((d2 - e) + (d1 - (t1 - e)));

	    /* The result is t1 + t2, after normalization. */
	    head_e1 = t1 + t2;
	    tail_e1 = t2 - (head_e1 - t1);
	  }
	  head_tmpx[0] = head_e1;
	  tail_tmpx[0] = tail_e1;
	  /* imaginary part */
	  d1 = (double) alpha_i[0] * x_ii[1];
	  d2 = (double) alpha_i[1] * x_ii[0];
	  {
	    /* Compute double-double = double + double. */
	    double e, t1, t2;

	    /* Knuth trick. */
	    t1 = d1 + d2;
	    e = t1 - d1;
	    t2 = ((d2 - e) + (d1 - (t1 - e)));

	    /* The result is t1 + t2, after normalization. */
	    head_e1 = t1 + t2;
	    tail_e1 = t2 - (head_e1 - t1);
	  }
	  head_tmpx[1] = head_e1;
	  tail_tmpx[1] = tail_e1;
	}			/* tmpx  = alpha * x[ix] */
	{
	  head_tmpy[0] = (double) beta_i[0] * y_ii;
	  tail_tmpy[0] = 0.0;
	  head_tmpy[1] = (double) beta_i[1] * y_ii;
	  tail_tmpy[1] = 0.0;
	}			/* tmpy = beta * y[iy] */
	{
	  double head_t, tail_t;
	  double head_a, tail_a;
	  double head_b, tail_b;
	  /* Real part */
	  head_a = head_tmpy[0];
	  tail_a = tail_tmpy[0];
	  head_b = head_tmpx[0];
	  tail_b = tail_tmpx[0];
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_a + head_b;
	    bv = s1 - head_a;
	    s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_a + tail_b;
	    bv = t1 - tail_a;
	    t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmpy[0] = head_t;
	  tail_tmpy[0] = tail_t;
	  /* Imaginary part */
	  head_a = head_tmpy[1];
	  tail_a = tail_tmpy[1];
	  head_b = head_tmpx[1];
	  tail_b = tail_tmpx[1];
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_a + head_b;
	    bv = s1 - head_a;
	    s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_a + tail_b;
	    bv = t1 - tail_a;
	    t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmpy[1] = head_t;
	  tail_tmpy[1] = tail_t;
	}
	w_i[iw] = head_tmpy[0];
	w_i[iw + 1] = head_tmpy[1];
	ix += incx;
	iy += incy;
	iw += incw;
      }				/* endfor */

      FPU_FIX_STOP;

      break;
    }
  }
}
