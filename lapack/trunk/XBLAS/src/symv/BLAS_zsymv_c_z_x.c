#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zsymv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a Symmetric matrix.
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
 * alpha   (input) const void*
 * 
 * a       (input) void*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x       (input) void*
 *         Vector x.
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) const void*
 * 
 * y       (input/output) void*
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
  static const char routine_name[] = "BLAS_zsymv_c_z_x";
  switch (prec) {

  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      /* Integer Index Variables */
      int i, k;

      int xi, yi;
      int aik, astarti, x_starti, y_starti;

      int incai;
      int incaik, incaik2;

      int n_i;

      /* Input Matrices */
      const float *a_i = (float *) a;
      const double *x_i = (double *) x;

      /* Output Vector */
      double *y_i = (double *) y;

      /* Input Scalars */
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;

      /* Temporary Floating-Point Variables */
      float a_elem[2];
      double x_elem[2];
      double y_elem[2];
      double prod[2];
      double sum[2];
      double tmp1[2];
      double tmp2[2];



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Check for error conditions. */
      if (lda < n) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -8, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -11, incy, NULL);
      }


      /* Set Index Parameters */
      n_i = n;

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaik = 1;
	incaik2 = lda;
      } else {
	incai = 1;
	incaik = lda;
	incaik2 = 1;
      }

      /* Adjustment to increments (if any) */
      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaik *= 2;
      incaik2 *= 2;
      if (incx < 0) {
	x_starti = (-n + 1) * incx;
      } else {
	x_starti = 0;
      }
      if (incy < 0) {
	y_starti = (-n + 1) * incy;
      } else {
	y_starti = 0;
      }



      /* alpha = 0.  In this case, just return beta * y */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp1[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp1[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  y_i[yi] = tmp1[0];
	  y_i[yi + 1] = tmp1[1];
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    sum[0] = sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti; k < i;
		 k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		prod[0] =
		  (double) a_elem[0] * x_elem[0] -
		  (double) a_elem[1] * x_elem[1];
		prod[1] =
		  (double) a_elem[0] * x_elem[1] +
		  (double) a_elem[1] * x_elem[0];
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		prod[0] =
		  (double) a_elem[0] * x_elem[0] -
		  (double) a_elem[1] * x_elem[1];
		prod[1] =
		  (double) a_elem[0] * x_elem[1] +
		  (double) a_elem[1] * x_elem[0];
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    y_i[yi] = sum[0];
	    y_i[yi + 1] = sum[1];
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    sum[0] = sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		prod[0] =
		  (double) a_elem[0] * x_elem[0] -
		  (double) a_elem[1] * x_elem[1];
		prod[1] =
		  (double) a_elem[0] * x_elem[1] +
		  (double) a_elem[1] * x_elem[0];
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		prod[0] =
		  (double) a_elem[0] * x_elem[0] -
		  (double) a_elem[1] * x_elem[1];
		prod[1] =
		  (double) a_elem[0] * x_elem[1] +
		  (double) a_elem[1] * x_elem[0];
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    y_elem[0] = y_i[yi];
	    y_elem[1] = y_i[yi + 1];
	    {
	      tmp2[0] =
		(double) y_elem[0] * beta_i[0] -
		(double) y_elem[1] * beta_i[1];
	      tmp2[1] =
		(double) y_elem[0] * beta_i[1] +
		(double) y_elem[1] * beta_i[0];
	    }
	    tmp1[0] = sum[0];
	    tmp1[1] = sum[1];
	    tmp1[0] = tmp2[0] + tmp1[0];
	    tmp1[1] = tmp2[1] + tmp1[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	  }
	}

      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti, astarti = 0;
	     i < n_i; i++, yi += incy, astarti += incai) {
	  sum[0] = sum[1] = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      prod[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      prod[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  {
	    tmp1[0] =
	      (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	  }
	  tmp1[0] = tmp2[0] + tmp1[0];
	  tmp1[1] = tmp2[1] + tmp1[1];
	  y_i[yi] = tmp1[0];
	  y_i[yi + 1] = tmp1[1];
	}
      }



      break;
    }

  case blas_prec_extra:{

      /* Integer Index Variables */
      int i, k;

      int xi, yi;
      int aik, astarti, x_starti, y_starti;

      int incai;
      int incaik, incaik2;

      int n_i;

      /* Input Matrices */
      const float *a_i = (float *) a;
      const double *x_i = (double *) x;

      /* Output Vector */
      double *y_i = (double *) y;

      /* Input Scalars */
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;

      /* Temporary Floating-Point Variables */
      float a_elem[2];
      double x_elem[2];
      double y_elem[2];
      double head_prod[2], tail_prod[2];
      double head_sum[2], tail_sum[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];

      FPU_FIX_DECL;

      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Check for error conditions. */
      if (lda < n) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -8, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -11, incy, NULL);
      }


      /* Set Index Parameters */
      n_i = n;

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaik = 1;
	incaik2 = lda;
      } else {
	incai = 1;
	incaik = lda;
	incaik2 = 1;
      }

      /* Adjustment to increments (if any) */
      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaik *= 2;
      incaik2 *= 2;
      if (incx < 0) {
	x_starti = (-n + 1) * incx;
      } else {
	x_starti = 0;
      }
      if (incy < 0) {
	y_starti = (-n + 1) * incy;
      } else {
	y_starti = 0;
      }

      FPU_FIX_START;

      /* alpha = 0.  In this case, just return beta * y */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[0] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[1] * beta_i[1];
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
	    head_tmp1[0] = head_t1;
	    tail_tmp1[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[1] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[0] * beta_i[1];
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
	    head_tmp1[1] = head_t1;
	    tail_tmp1[1] = tail_t1;
	  }
	  y_i[yi] = head_tmp1[0];
	  y_i[yi + 1] = head_tmp1[1];
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti; k < i;
		 k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		double cd[2];
		cd[0] = (double) a_elem[0];
		cd[1] = (double) a_elem[1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[0] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[1] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[0] = head_t1;
		  tail_prod[0] = tail_t1;
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[1] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[0] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[1] = head_t1;
		  tail_prod[1] = tail_t1;
		}
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
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
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
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
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		double cd[2];
		cd[0] = (double) a_elem[0];
		cd[1] = (double) a_elem[1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[0] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[1] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[0] = head_t1;
		  tail_prod[0] = tail_t1;
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[1] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[0] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[1] = head_t1;
		  tail_prod[1] = tail_t1;
		}
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
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
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
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
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }
	    }
	    y_i[yi] = head_sum[0];
	    y_i[yi + 1] = head_sum[1];
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		double cd[2];
		cd[0] = (double) a_elem[0];
		cd[1] = (double) a_elem[1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[0] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[1] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[0] = head_t1;
		  tail_prod[0] = tail_t1;
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[1] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[0] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[1] = head_t1;
		  tail_prod[1] = tail_t1;
		}
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
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
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
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
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		double cd[2];
		cd[0] = (double) a_elem[0];
		cd[1] = (double) a_elem[1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[0] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[1] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[0] = head_t1;
		  tail_prod[0] = tail_t1;
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[1] * split;
		    a1 = con - x_elem[1];
		    a1 = con - a1;
		    a2 = x_elem[1] - a1;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    head_t1 = x_elem[1] * cd[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = x_elem[0] * split;
		    a1 = con - x_elem[0];
		    a1 = con - a1;
		    a2 = x_elem[0] - a1;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    head_t2 = x_elem[0] * cd[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  head_prod[1] = head_t1;
		  tail_prod[1] = tail_t1;
		}
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
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
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
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
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }
	    }
	    y_elem[0] = y_i[yi];
	    y_elem[1] = y_i[yi + 1];
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[0] * split;
		a1 = con - y_elem[0];
		a1 = con - a1;
		a2 = y_elem[0] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = y_elem[0] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[1] * split;
		a1 = con - y_elem[1];
		a1 = con - a1;
		a2 = y_elem[1] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = y_elem[1] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      head_tmp2[0] = head_t1;
	      tail_tmp2[0] = tail_t1;
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[1] * split;
		a1 = con - y_elem[1];
		a1 = con - a1;
		a2 = y_elem[1] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = y_elem[1] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[0] * split;
		a1 = con - y_elem[0];
		a1 = con - a1;
		a2 = y_elem[0] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = y_elem[0] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      head_tmp2[1] = head_t1;
	      tail_tmp2[1] = tail_t1;
	    }
	    head_tmp1[0] = head_sum[0];
	    tail_tmp1[0] = tail_sum[0];
	    head_tmp1[1] = head_sum[1];
	    tail_tmp1[1] = tail_sum[1];
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp2[0];
	      tail_a = tail_tmp2[0];
	      head_b = head_tmp1[0];
	      tail_b = tail_tmp1[0];
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
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp2[1];
	      tail_a = tail_tmp2[1];
	      head_b = head_tmp1[1];
	      tail_b = tail_tmp1[1];
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
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	  }
	}

      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti, astarti = 0;
	     i < n_i; i++, yi += incy, astarti += incai) {
	  head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      double cd[2];
	      cd[0] = (double) a_elem[0];
	      cd[1] = (double) a_elem[1];
	      {
		/* Compute complex-extra = complex-double * complex-double. */
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		/* Real part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[0] * split;
		  a1 = con - x_elem[0];
		  a1 = con - a1;
		  a2 = x_elem[0] - a1;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  head_t1 = x_elem[0] * cd[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[1] * split;
		  a1 = con - x_elem[1];
		  a1 = con - a1;
		  a2 = x_elem[1] - a1;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  head_t2 = x_elem[1] * cd[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_prod[0] = head_t1;
		tail_prod[0] = tail_t1;
		/* Imaginary part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[1] * split;
		  a1 = con - x_elem[1];
		  a1 = con - a1;
		  a2 = x_elem[1] - a1;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  head_t1 = x_elem[1] * cd[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[0] * split;
		  a1 = con - x_elem[0];
		  a1 = con - a1;
		  a2 = x_elem[0] - a1;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  head_t2 = x_elem[0] * cd[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_prod[1] = head_t1;
		tail_prod[1] = tail_t1;
	      }
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum[0];
	      tail_a = tail_sum[0];
	      head_b = head_prod[0];
	      tail_b = tail_prod[0];
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
	      head_sum[0] = head_t;
	      tail_sum[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum[1];
	      tail_a = tail_sum[1];
	      head_b = head_prod[1];
	      tail_b = tail_prod[1];
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
	      head_sum[1] = head_t;
	      tail_sum[1] = tail_t;
	    }
	  }
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      double cd[2];
	      cd[0] = (double) a_elem[0];
	      cd[1] = (double) a_elem[1];
	      {
		/* Compute complex-extra = complex-double * complex-double. */
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		/* Real part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[0] * split;
		  a1 = con - x_elem[0];
		  a1 = con - a1;
		  a2 = x_elem[0] - a1;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  head_t1 = x_elem[0] * cd[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[1] * split;
		  a1 = con - x_elem[1];
		  a1 = con - a1;
		  a2 = x_elem[1] - a1;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  head_t2 = x_elem[1] * cd[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_prod[0] = head_t1;
		tail_prod[0] = tail_t1;
		/* Imaginary part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[1] * split;
		  a1 = con - x_elem[1];
		  a1 = con - a1;
		  a2 = x_elem[1] - a1;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  head_t1 = x_elem[1] * cd[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem[0] * split;
		  a1 = con - x_elem[0];
		  a1 = con - a1;
		  a2 = x_elem[0] - a1;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  head_t2 = x_elem[0] * cd[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_prod[1] = head_t1;
		tail_prod[1] = tail_t1;
	      }
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum[0];
	      tail_a = tail_sum[0];
	      head_b = head_prod[0];
	      tail_b = tail_prod[0];
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
	      head_sum[0] = head_t;
	      tail_sum[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum[1];
	      tail_a = tail_sum[1];
	      head_b = head_prod[1];
	      tail_b = tail_prod[1];
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
	      head_sum[1] = head_t;
	      tail_sum[1] = tail_t;
	    }
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[0] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[1] * beta_i[1];
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
	    head_tmp2[0] = head_t1;
	    tail_tmp2[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[1] * split;
	      a1 = con - y_elem[1];
	      a1 = con - a1;
	      a2 = y_elem[1] - a1;
	      con = beta_i[0] * split;
	      b1 = con - beta_i[0];
	      b1 = con - b1;
	      b2 = beta_i[0] - b1;

	      head_t1 = y_elem[1] * beta_i[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = y_elem[0] * split;
	      a1 = con - y_elem[0];
	      a1 = con - a1;
	      a2 = y_elem[0] - a1;
	      con = beta_i[1] * split;
	      b1 = con - beta_i[1];
	      b1 = con - b1;
	      b2 = beta_i[1] - b1;

	      head_t2 = y_elem[0] * beta_i[1];
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
	    head_tmp2[1] = head_t1;
	    tail_tmp2[1] = tail_t1;
	  }
	  {
	    /* Compute complex-extra = complex-extra * complex-double. */
	    double head_a0, tail_a0;
	    double head_a1, tail_a1;
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    head_a0 = head_sum[0];
	    tail_a0 = tail_sum[0];
	    head_a1 = head_sum[1];
	    tail_a1 = tail_sum[1];
	    /* real part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a0 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a1 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
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
	    head_tmp1[0] = head_t1;
	    tail_tmp1[0] = tail_t1;
	    /* imaginary part */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a1 * split;
	      a11 = con - head_a1;
	      a11 = con - a11;
	      a21 = head_a1 - a11;
	      con = alpha_i[0] * split;
	      b1 = con - alpha_i[0];
	      b1 = con - b1;
	      b2 = alpha_i[0] - b1;

	      c11 = head_a1 * alpha_i[0];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a1 * alpha_i[0];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_a0 * split;
	      a11 = con - head_a0;
	      a11 = con - a11;
	      a21 = head_a0 - a11;
	      con = alpha_i[1] * split;
	      b1 = con - alpha_i[1];
	      b1 = con - b1;
	      b2 = alpha_i[1] - b1;

	      c11 = head_a0 * alpha_i[1];
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_a0 * alpha_i[1];
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_t2 = t1 + t2;
	      tail_t2 = t2 - (head_t2 - t1);
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
	    head_tmp1[1] = head_t1;
	    tail_tmp1[1] = tail_t1;
	  }

	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_tmp2[0];
	    tail_a = tail_tmp2[0];
	    head_b = head_tmp1[0];
	    tail_b = tail_tmp1[0];
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
	    head_tmp1[0] = head_t;
	    tail_tmp1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_tmp2[1];
	    tail_a = tail_tmp2[1];
	    head_b = head_tmp1[1];
	    tail_b = tail_tmp1[1];
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
	    head_tmp1[1] = head_t;
	    tail_tmp1[1] = tail_t;
	  }
	  y_i[yi] = head_tmp1[0];
	  y_i[yi + 1] = head_tmp1[1];
	}
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_zsymv_c_z_x */
