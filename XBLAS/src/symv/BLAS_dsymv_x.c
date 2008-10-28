#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dsymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, double alpha, const double *a, int lda,
		  const double *x, int incx, double beta,
		  double *y, int incy, enum blas_prec_type prec)

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
 * alpha   (input) double
 * 
 * a       (input) double*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x       (input) double*
 *         Vector x.
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) double
 * 
 * y       (input/output) double*
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
  static const char routine_name[] = "BLAS_dsymv_x";
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
      const double *a_i = a;
      const double *x_i = x;

      /* Output Vector */
      double *y_i = y;

      /* Input Scalars */
      double alpha_i = alpha;
      double beta_i = beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
      double x_elem;
      double y_elem;
      double prod;
      double sum;
      double tmp1;
      double tmp2;



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
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
      if (alpha_i == 0.0) {
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  y_elem = y_i[yi];
	  tmp1 = y_elem * beta_i;
	  y_i[yi] = tmp1;
	}
      } else if (alpha_i == 1.0) {

	/* Case alpha == 1. */

	if (beta_i == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    sum = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti; k < i;
		 k++, aik += incaik, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	    }
	    y_i[yi] = sum;
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    sum = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	    }
	    y_elem = y_i[yi];
	    tmp2 = y_elem * beta_i;
	    tmp1 = sum;
	    tmp1 = tmp2 + tmp1;
	    y_i[yi] = tmp1;
	  }
	}

      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti, astarti = 0;
	     i < n_i; i++, yi += incy, astarti += incai) {
	  sum = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem = a_i[aik];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem = a_i[aik];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  y_elem = y_i[yi];
	  tmp2 = y_elem * beta_i;
	  tmp1 = sum * alpha_i;
	  tmp1 = tmp2 + tmp1;
	  y_i[yi] = tmp1;
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
      const double *a_i = a;
      const double *x_i = x;

      /* Output Vector */
      double *y_i = y;

      /* Input Scalars */
      double alpha_i = alpha;
      double beta_i = beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
      double x_elem;
      double y_elem;
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;

      FPU_FIX_DECL;

      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
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
      if (alpha_i == 0.0) {
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
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

	    head_tmp1 = y_elem * beta_i;
	    tail_tmp1 =
	      (((a1 * b1 - head_tmp1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  y_i[yi] = head_tmp1;
	}
      } else if (alpha_i == 1.0) {

	/* Case alpha == 1. */

	if (beta_i == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    head_sum = tail_sum = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti; k < i;
		 k++, aik += incaik, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = a_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = a_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	    }
	    y_i[yi] = head_sum;
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incai) {
	    head_sum = tail_sum = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = a_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	    }
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem = a_i[aik];
	      x_elem = x_i[xi];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = a_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
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
	      tail_tmp2 =
		(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_tmp1 = head_sum;
	    tail_tmp1 = tail_sum;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_tmp2 + head_tmp1;
	      bv = s1 - head_tmp2;
	      s2 = ((head_tmp1 - bv) + (head_tmp2 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_tmp2 + tail_tmp1;
	      bv = t1 - tail_tmp2;
	      t2 = ((tail_tmp1 - bv) + (tail_tmp2 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	  }
	}

      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti, astarti = 0;
	     i < n_i; i++, yi += incy, astarti += incai) {
	  head_sum = tail_sum = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem = a_i[aik];
	    x_elem = x_i[xi];
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = a_elem * split;
	      a1 = con - a_elem;
	      a1 = con - a1;
	      a2 = a_elem - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod = a_elem * x_elem;
	      tail_prod =
		(((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_sum + head_prod;
	      bv = s1 - head_sum;
	      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_sum + tail_prod;
	      bv = t1 - tail_sum;
	      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_sum = t1 + t2;
	      tail_sum = t2 - (head_sum - t1);
	    }
	  }
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem = a_i[aik];
	    x_elem = x_i[xi];
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = a_elem * split;
	      a1 = con - a_elem;
	      a1 = con - a1;
	      a2 = a_elem - a1;
	      con = x_elem * split;
	      b1 = con - x_elem;
	      b1 = con - b1;
	      b2 = x_elem - b1;

	      head_prod = a_elem * x_elem;
	      tail_prod =
		(((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_sum + head_prod;
	      bv = s1 - head_sum;
	      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_sum + tail_prod;
	      bv = t1 - tail_sum;
	      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_sum = t1 + t2;
	      tail_sum = t2 - (head_sum - t1);
	    }
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
	    tail_tmp2 =
	      (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_sum * split;
	    a11 = con - head_sum;
	    a11 = con - a11;
	    a21 = head_sum - a11;
	    con = alpha_i * split;
	    b1 = con - alpha_i;
	    b1 = con - b1;
	    b2 = alpha_i - b1;

	    c11 = head_sum * alpha_i;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_sum * alpha_i;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_tmp1 = t1 + t2;
	    tail_tmp1 = t2 - (head_tmp1 - t1);
	  }
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_tmp2 + head_tmp1;
	    bv = s1 - head_tmp2;
	    s2 = ((head_tmp1 - bv) + (head_tmp2 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_tmp2 + tail_tmp1;
	    bv = t1 - tail_tmp2;
	    t2 = ((tail_tmp1 - bv) + (tail_tmp2 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_tmp1 = t1 + t2;
	    tail_tmp1 = t2 - (head_tmp1 - t1);
	  }
	  y_i[yi] = head_tmp1;
	}
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_dsymv_x */
