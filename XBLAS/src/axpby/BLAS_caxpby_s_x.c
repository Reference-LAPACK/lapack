#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_caxpby_s_x(int n, const void *alpha, const float *x, int incx,
		     const void *beta, void *y,
		     int incy, enum blas_prec_type prec)

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
 * x         (input) const float*
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
  static const char routine_name[] = "BLAS_caxpby_s_x";

  switch (prec) {
  case blas_prec_single:{

      int i, ix = 0, iy = 0;
      const float *x_i = x;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii[2];
      float tmpx[2];
      float tmpy[2];


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




      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii[0] = y_i[iy];
	y_ii[1] = y_i[iy + 1];
	{
	  tmpx[0] = alpha_i[0] * x_ii;
	  tmpx[1] = alpha_i[1] * x_ii;
	}			/* tmpx  = alpha * x[ix] */
	{
	  tmpy[0] = beta_i[0] * y_ii[0] - beta_i[1] * y_ii[1];
	  tmpy[1] = beta_i[0] * y_ii[1] + beta_i[1] * y_ii[0];
	}
	/* tmpy = beta * y[iy] */
	tmpy[0] = tmpy[0] + tmpx[0];
	tmpy[1] = tmpy[1] + tmpx[1];
	y_i[iy] = tmpy[0];
	y_i[iy + 1] = tmpy[1];
	ix += incx;
	iy += incy;
      }				/* endfor */



      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      const float *x_i = x;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii[2];
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




      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii[0] = y_i[iy];
	y_ii[1] = y_i[iy + 1];
	{
	  tmpx[0] = (double) alpha_i[0] * x_ii;
	  tmpx[1] = (double) alpha_i[1] * x_ii;
	}			/* tmpx  = alpha * x[ix] */
	{
	  tmpy[0] =
	    (double) beta_i[0] * y_ii[0] - (double) beta_i[1] * y_ii[1];
	  tmpy[1] =
	    (double) beta_i[0] * y_ii[1] + (double) beta_i[1] * y_ii[0];
	}			/* tmpy = beta * y[iy] */
	tmpy[0] = tmpy[0] + tmpx[0];
	tmpy[1] = tmpy[1] + tmpx[1];
	y_i[iy] = tmpy[0];
	y_i[iy + 1] = tmpy[1];
	ix += incx;
	iy += incy;
      }				/* endfor */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      const float *x_i = x;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii[2];
      double head_tmpx[2], tail_tmpx[2];
      double head_tmpy[2], tail_tmpy[2];
      FPU_FIX_DECL;

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

      FPU_FIX_START;


      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii[0] = y_i[iy];
	y_ii[1] = y_i[iy + 1];
	{
	  head_tmpx[0] = (double) alpha_i[0] * x_ii;
	  tail_tmpx[0] = 0.0;
	  head_tmpx[1] = (double) alpha_i[1] * x_ii;
	  tail_tmpx[1] = 0.0;
	}			/* tmpx  = alpha * x[ix] */
	{
	  double head_e1, tail_e1;
	  double d1;
	  double d2;
	  /* Real part */
	  d1 = (double) beta_i[0] * y_ii[0];
	  d2 = (double) -beta_i[1] * y_ii[1];
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
	  head_tmpy[0] = head_e1;
	  tail_tmpy[0] = tail_e1;
	  /* imaginary part */
	  d1 = (double) beta_i[0] * y_ii[1];
	  d2 = (double) beta_i[1] * y_ii[0];
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
	  head_tmpy[1] = head_e1;
	  tail_tmpy[1] = tail_e1;
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
	y_i[iy] = head_tmpy[0];
	y_i[iy + 1] = head_tmpy[1];
	ix += incx;
	iy += incy;
      }				/* endfor */

      FPU_FIX_STOP;
    }
    break;
  }
}				/* end BLAS_caxpby_s_x */
