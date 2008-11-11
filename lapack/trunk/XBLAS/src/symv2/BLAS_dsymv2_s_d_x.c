#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_dsymv2_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double alpha, const float *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a symmetric matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input symmetric matrix A.
 * 
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * n       (input) int
 *         Dimension of A and size of vectors x, y.
 *
 * alpha   (input) double
 * 
 * a       (input) float*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) double*
 *         Vector x_head
 *
 * x_tail  (input) double*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) double
 * 
 * y       (input) double*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
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
  /* Routine name */
  const char routine_name[] = "BLAS_dsymv2_s_d_x";
  switch (prec) {

  case blas_prec_single:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const float *a_i = a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float a_elem;
      double x_elem;
      double y_elem;
      double prod1;
      double prod2;
      double sum1;
      double sum2;
      double tmp1;
      double tmp2;
      double tmp3;



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	sum1 = 0.0;
	sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	sum1 = sum1 + sum2;
	tmp1 = sum1 * alpha_i;
	y_elem = y_i[yi];
	tmp2 = y_elem * beta_i;
	tmp3 = tmp1 + tmp2;
	y_i[yi] = tmp3;
      }



      break;
    }

  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const float *a_i = a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float a_elem;
      double x_elem;
      double y_elem;
      double prod1;
      double prod2;
      double sum1;
      double sum2;
      double tmp1;
      double tmp2;
      double tmp3;



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	sum1 = 0.0;
	sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	sum1 = sum1 + sum2;
	tmp1 = sum1 * alpha_i;
	y_elem = y_i[yi];
	tmp2 = y_elem * beta_i;
	tmp3 = tmp1 + tmp2;
	y_i[yi] = tmp3;
      }



      break;
    }

  case blas_prec_extra:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const float *a_i = a;
      const double *x_head_i = x_head;
      const double *x_tail_i = x_tail;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float a_elem;
      double x_elem;
      double y_elem;
      double head_prod1, tail_prod1;
      double head_prod2, tail_prod2;
      double head_sum1, tail_sum1;
      double head_sum2, tail_sum2;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;
      double head_tmp3, tail_tmp3;

      FPU_FIX_DECL;

      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);

      FPU_FIX_START;

      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	head_sum1 = tail_sum1 = 0.0;
	head_sum2 = tail_sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    double dt = (double) a_elem;
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = dt * split;
	      a1 = con - dt;
	      a1 = con - a1;
	      a2 = dt - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod1 = dt * x_elem;
	      tail_prod1 =
		(((a1 * b1 - head_prod1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum1 + head_prod1;
	    bv = s1 - head_sum1;
	    s2 = ((head_prod1 - bv) + (head_sum1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum1 + tail_prod1;
	    bv = t1 - tail_sum1;
	    t2 = ((tail_prod1 - bv) + (tail_sum1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum1 = t1 + t2;
	    tail_sum1 = t2 - (head_sum1 - t1);
	  }
	  x_elem = x_tail_i[xi];
	  {
	    double dt = (double) a_elem;
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = dt * split;
	      a1 = con - dt;
	      a1 = con - a1;
	      a2 = dt - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod2 = dt * x_elem;
	      tail_prod2 =
		(((a1 * b1 - head_prod2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum2 + head_prod2;
	    bv = s1 - head_sum2;
	    s2 = ((head_prod2 - bv) + (head_sum2 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum2 + tail_prod2;
	    bv = t1 - tail_sum2;
	    t2 = ((tail_prod2 - bv) + (tail_sum2 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum2 = t1 + t2;
	    tail_sum2 = t2 - (head_sum2 - t1);
	  }
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    double dt = (double) a_elem;
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = dt * split;
	      a1 = con - dt;
	      a1 = con - a1;
	      a2 = dt - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod1 = dt * x_elem;
	      tail_prod1 =
		(((a1 * b1 - head_prod1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum1 + head_prod1;
	    bv = s1 - head_sum1;
	    s2 = ((head_prod1 - bv) + (head_sum1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum1 + tail_prod1;
	    bv = t1 - tail_sum1;
	    t2 = ((tail_prod1 - bv) + (tail_sum1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum1 = t1 + t2;
	    tail_sum1 = t2 - (head_sum1 - t1);
	  }
	  x_elem = x_tail_i[xi];
	  {
	    double dt = (double) a_elem;
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = dt * split;
	      a1 = con - dt;
	      a1 = con - a1;
	      a2 = dt - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod2 = dt * x_elem;
	      tail_prod2 =
		(((a1 * b1 - head_prod2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum2 + head_prod2;
	    bv = s1 - head_sum2;
	    s2 = ((head_prod2 - bv) + (head_sum2 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum2 + tail_prod2;
	    bv = t1 - tail_sum2;
	    t2 = ((tail_prod2 - bv) + (tail_sum2 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum2 = t1 + t2;
	    tail_sum2 = t2 - (head_sum2 - t1);
	  }
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_sum1 + head_sum2;
	  bv = s1 - head_sum1;
	  s2 = ((head_sum2 - bv) + (head_sum1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_sum1 + tail_sum2;
	  bv = t1 - tail_sum1;
	  t2 = ((tail_sum2 - bv) + (tail_sum1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_sum1 = t1 + t2;
	  tail_sum1 = t2 - (head_sum1 - t1);
	}
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_sum1 * split;
	  a11 = con - head_sum1;
	  a11 = con - a11;
	  a21 = head_sum1 - a11;
	  con = alpha_i * split;
	  b1 = con - alpha_i;
	  b1 = con - b1;
	  b2 = alpha_i - b1;

	  c11 = head_sum1 * alpha_i;
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_sum1 * alpha_i;
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_tmp1 = t1 + t2;
	  tail_tmp1 = t2 - (head_tmp1 - t1);
	}
	y_elem = y_i[yi];
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = y_elem * split;
	  a1 = con - y_elem;
	  a1 = con - a1;
	  a2 = y_elem - a1;
	  con = beta_i * split;
	  b1 = con - beta_i;
	  b1 = con - b1;
	  b2 = beta_i - b1;

	  head_tmp2 = y_elem * beta_i;
	  tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_tmp1 + head_tmp2;
	  bv = s1 - head_tmp1;
	  s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_tmp1 + tail_tmp2;
	  bv = t1 - tail_tmp1;
	  t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_tmp3 = t1 + t2;
	  tail_tmp3 = t2 - (head_tmp3 - t1);
	}
	y_i[yi] = head_tmp3;
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_dsymv2_s_d_x */
