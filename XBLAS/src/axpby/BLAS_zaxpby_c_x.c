#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zaxpby_c_x(int n, const void *alpha, const void *x, int incx,
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
  static const char routine_name[] = "BLAS_zaxpby_c_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
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
	  tmpx[0] =
	    (double) alpha_i[0] * x_ii[0] - (double) alpha_i[1] * x_ii[1];
	  tmpx[1] =
	    (double) alpha_i[0] * x_ii[1] + (double) alpha_i[1] * x_ii[0];
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
      const float *x_i = (float *) x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      float x_ii[2];
      double y_ii[2];
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
	  double cd[2];
	  cd[0] = (double) x_ii[0];
	  cd[1] = (double) x_ii[1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = alpha_i[0] * split;
	      a1 = con - alpha_i[0];
	      a1 = con - a1;
	      a2 = alpha_i[0] - a1;
	      con = cd[0] * split;
	      b1 = con - cd[0];
	      b1 = con - b1;
	      b2 = cd[0] - b1;

	      head_t1 = alpha_i[0] * cd[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = alpha_i[1] * split;
	      a1 = con - alpha_i[1];
	      a1 = con - a1;
	      a2 = alpha_i[1] - a1;
	      con = cd[1] * split;
	      b1 = con - cd[1];
	      b1 = con - b1;
	      b2 = cd[1] - b1;

	      head_t2 = alpha_i[1] * cd[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmpx[0] = head_t1;
	    tail_tmpx[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = alpha_i[1] * split;
	      a1 = con - alpha_i[1];
	      a1 = con - a1;
	      a2 = alpha_i[1] - a1;
	      con = cd[0] * split;
	      b1 = con - cd[0];
	      b1 = con - b1;
	      b2 = cd[0] - b1;

	      head_t1 = alpha_i[1] * cd[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = alpha_i[0] * split;
	      a1 = con - alpha_i[0];
	      a1 = con - a1;
	      a2 = alpha_i[0] - a1;
	      con = cd[1] * split;
	      b1 = con - cd[1];
	      b1 = con - b1;
	      b2 = cd[1] - b1;

	      head_t2 = alpha_i[0] * cd[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_tmpx[1] = head_t1;
	    tail_tmpx[1] = tail_t1;
	  }
	}			/* tmpx  = alpha * x[ix] */
	{
	  /* Compute complex-extra = complex-double * complex-double. */
	  double head_t1, tail_t1;
	  double head_t2, tail_t2;
	  /* Real part */
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_ii[0] * split;
	    b1 = con - y_ii[0];
	    b1 = con - b1;
	    b2 = y_ii[0] - b1;

	    head_t1 = beta_i[0] * y_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_ii[1] * split;
	    b1 = con - y_ii[1];
	    b1 = con - b1;
	    b2 = y_ii[1] - b1;

	    head_t2 = beta_i[1] * y_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  head_t2 = -head_t2;
	  tail_t2 = -tail_t2;
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpy[0] = head_t1;
	  tail_tmpy[0] = tail_t1;
	  /* Imaginary part */
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_ii[0] * split;
	    b1 = con - y_ii[0];
	    b1 = con - b1;
	    b2 = y_ii[0] - b1;

	    head_t1 = beta_i[1] * y_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_ii[1] * split;
	    b1 = con - y_ii[1];
	    b1 = con - b1;
	    b2 = y_ii[1] - b1;

	    head_t2 = beta_i[0] * y_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpy[1] = head_t1;
	  tail_tmpy[1] = tail_t1;
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
}				/* end BLAS_zaxpby_c_x */
