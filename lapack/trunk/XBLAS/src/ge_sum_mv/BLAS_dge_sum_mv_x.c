#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dge_sum_mv_x(enum blas_order_type order, int m, int n,
		       double alpha, const double *a, int lda,
		       const double *x, int incx,
		       double beta, const double *b, int ldb,
		       double *y, int incy, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * B * y, 
 *     where A, B are general matricies.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * m            (input) int
 *              Row Dimension of A, B, length of output vector y
 *
 * n            (input) int
 *              Column Dimension of A, B and the length of vector x
 *
 * alpha        (input) double
 *              
 * A            (input) const double*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) double
 *
 * b            (input) const double*
 *
 * ldb          (input) int 
 *              Leading dimension of B
 *
 * y            (input/output) double*
 *
 * incy         (input) int
 *              The stride for vector y.
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
  static const char routine_name[] = "BLAS_dge_sum_mv";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_indigenous:
  case blas_prec_double:
    {
      int i, j;
      int xi, yi;
      int x_starti, y_starti, incxi, incyi;
      int lda_min;
      int ai;
      int incai;
      int aij;
      int incaij;
      int bi;
      int incbi;
      int bij;
      int incbij;

      const double *a_i = a;
      const double *b_i = b;
      const double *x_i = x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      double a_elem;
      double b_elem;
      double x_elem;
      double prod;
      double sumA;
      double sumB;
      double tmp1;
      double tmp2;



      /* m is number of rows */
      /* n is number of columns */

      if (m == 0 || n == 0)
	return;


      /* all error calls */
      if (order == blas_rowmajor) {
	lda_min = n;
	incai = lda;		/* row stride */
	incbi = ldb;
	incbij = incaij = 1;	/* column stride */
      } else if (order == blas_colmajor) {
	lda_min = m;
	incai = incbi = 1;	/*row stride */
	incaij = lda;		/* column stride */
	incbij = ldb;
      } else {
	/* error, order not blas_colmajor not blas_rowmajor */
	BLAS_error(routine_name, -1, order, 0);
	return;
      }

      if (m < 0)
	BLAS_error(routine_name, -2, m, 0);
      else if (n < 0)
	BLAS_error(routine_name, -3, n, 0);
      if (lda < lda_min)
	BLAS_error(routine_name, -6, lda, 0);
      else if (ldb < lda_min)
	BLAS_error(routine_name, -11, ldb, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -8, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      incxi = incx;
      incyi = incy;







      if (incxi > 0)
	x_starti = 0;
      else
	x_starti = (1 - n) * incxi;

      if (incyi > 0)
	y_starti = 0;
      else
	y_starti = (1 - m) * incyi;



      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  /* alpha, beta are 0.0 */
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    y_i[yi] = 0.0;
	  }
	} else if (beta_i == 1.0) {
	  /* alpha is 0.0, beta is 1.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    y_i[yi] = sumB;

	    bi += incbi;
	  }
	} else {
	  /* alpha is 0.0, beta not 1.0 nor 0.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1 = sumB * beta_i;
	    y_i[yi] = tmp1;

	    bi += incbi;
	  }
	}
      } else if (alpha_i == 1.0) {
	if (beta_i == 0.0) {
	  /* alpha is 1.0, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    y_i[yi] = sumA;
	    ai += incai;

	  }
	} else if (beta_i == 1.0) {
	  /* alpha is 1.0, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;
	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;
	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1 = sumA;
	    tmp2 = sumB;
	    tmp1 = tmp1 + tmp2;
	    y_i[yi] = tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* alpha is 1.0, beta is other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;
	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;
	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1 = sumA;
	    tmp2 = sumB * beta_i;
	    tmp1 = tmp1 + tmp2;
	    y_i[yi] = tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	}
      } else {
	if (beta_i == 0.0) {
	  /* alpha is other, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    tmp1 = sumA * alpha_i;
	    y_i[yi] = tmp1;
	    ai += incai;

	  }
	} else if (beta_i == 1.0) {
	  /* alpha is other, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;
	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;
	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1 = sumA * alpha_i;
	    tmp2 = sumB;
	    tmp1 = tmp1 + tmp2;
	    y_i[yi] = tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* most general form, alpha, beta are other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA = 0.0;
	    aij = ai;
	    sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
	      prod = a_elem * x_elem;
	      sumA = sumA + prod;
	      aij += incaij;
	      b_elem = b_i[bij];
	      prod = b_elem * x_elem;
	      sumB = sumB + prod;
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1 = sumA * alpha_i;
	    tmp2 = sumB * beta_i;
	    tmp1 = tmp1 + tmp2;
	    y_i[yi] = tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	}
      }

    }
    break;

  case blas_prec_extra:
    {
      int i, j;
      int xi, yi;
      int x_starti, y_starti, incxi, incyi;
      int lda_min;
      int ai;
      int incai;
      int aij;
      int incaij;
      int bi;
      int incbi;
      int bij;
      int incbij;

      const double *a_i = a;
      const double *b_i = b;
      const double *x_i = x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      double a_elem;
      double b_elem;
      double x_elem;
      double head_prod, tail_prod;
      double head_sumA, tail_sumA;
      double head_sumB, tail_sumB;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;

      FPU_FIX_DECL;

      /* m is number of rows */
      /* n is number of columns */

      if (m == 0 || n == 0)
	return;


      /* all error calls */
      if (order == blas_rowmajor) {
	lda_min = n;
	incai = lda;		/* row stride */
	incbi = ldb;
	incbij = incaij = 1;	/* column stride */
      } else if (order == blas_colmajor) {
	lda_min = m;
	incai = incbi = 1;	/*row stride */
	incaij = lda;		/* column stride */
	incbij = ldb;
      } else {
	/* error, order not blas_colmajor not blas_rowmajor */
	BLAS_error(routine_name, -1, order, 0);
	return;
      }

      if (m < 0)
	BLAS_error(routine_name, -2, m, 0);
      else if (n < 0)
	BLAS_error(routine_name, -3, n, 0);
      if (lda < lda_min)
	BLAS_error(routine_name, -6, lda, 0);
      else if (ldb < lda_min)
	BLAS_error(routine_name, -11, ldb, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -8, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      incxi = incx;
      incyi = incy;







      if (incxi > 0)
	x_starti = 0;
      else
	x_starti = (1 - n) * incxi;

      if (incyi > 0)
	y_starti = 0;
      else
	y_starti = (1 - m) * incyi;

      FPU_FIX_START;

      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  /* alpha, beta are 0.0 */
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    y_i[yi] = 0.0;
	  }
	} else if (beta_i == 1.0) {
	  /* alpha is 0.0, beta is 1.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    y_i[yi] = head_sumB;

	    bi += incbi;
	  }
	} else {
	  /* alpha is 0.0, beta not 1.0 nor 0.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumB * split;
	      a11 = con - head_sumB;
	      a11 = con - a11;
	      a21 = head_sumB - a11;
	      con = beta_i * split;
	      b1 = con - beta_i;
	      b1 = con - b1;
	      b2 = beta_i - b1;

	      c11 = head_sumB * beta_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumB * beta_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;

	    bi += incbi;
	  }
	}
      } else if (alpha_i == 1.0) {
	if (beta_i == 0.0) {
	  /* alpha is 1.0, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    y_i[yi] = head_sumA;
	    ai += incai;

	  }
	} else if (beta_i == 1.0) {
	  /* alpha is 1.0, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;
	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;
	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    head_tmp1 = head_sumA;
	    tail_tmp1 = tail_sumA;
	    head_tmp2 = head_sumB;
	    tail_tmp2 = tail_sumB;
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
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* alpha is 1.0, beta is other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;
	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;
	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    head_tmp1 = head_sumA;
	    tail_tmp1 = tail_sumA;
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumB * split;
	      a11 = con - head_sumB;
	      a11 = con - a11;
	      a21 = head_sumB - a11;
	      con = beta_i * split;
	      b1 = con - beta_i;
	      b1 = con - b1;
	      b2 = beta_i - b1;

	      c11 = head_sumB * beta_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumB * beta_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp2 = t1 + t2;
	      tail_tmp2 = t2 - (head_tmp2 - t1);
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
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	}
      } else {
	if (beta_i == 0.0) {
	  /* alpha is other, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumA * split;
	      a11 = con - head_sumA;
	      a11 = con - a11;
	      a21 = head_sumA - a11;
	      con = alpha_i * split;
	      b1 = con - alpha_i;
	      b1 = con - b1;
	      b2 = alpha_i - b1;

	      c11 = head_sumA * alpha_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumA * alpha_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	    ai += incai;

	  }
	} else if (beta_i == 1.0) {
	  /* alpha is other, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;
	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;
	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumA * split;
	      a11 = con - head_sumA;
	      a11 = con - a11;
	      a21 = head_sumA - a11;
	      con = alpha_i * split;
	      b1 = con - alpha_i;
	      b1 = con - b1;
	      b2 = alpha_i - b1;

	      c11 = head_sumA * alpha_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumA * alpha_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    head_tmp2 = head_sumB;
	    tail_tmp2 = tail_sumB;
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
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* most general form, alpha, beta are other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA = tail_sumA = 0.0;
	    aij = ai;
	    head_sumB = tail_sumB = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem = a_i[aij];
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
		s1 = head_sumA + head_prod;
		bv = s1 - head_sumA;
		s2 = ((head_prod - bv) + (head_sumA - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumA + tail_prod;
		bv = t1 - tail_sumA;
		t2 = ((tail_prod - bv) + (tail_sumA - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumA = t1 + t2;
		tail_sumA = t2 - (head_sumA - t1);
	      }
	      aij += incaij;
	      b_elem = b_i[bij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = b_elem * split;
		a1 = con - b_elem;
		a1 = con - a1;
		a2 = b_elem - a1;
		con = x_elem * split;
		b1 = con - x_elem;
		b1 = con - b1;
		b2 = x_elem - b1;

		head_prod = b_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sumB + head_prod;
		bv = s1 - head_sumB;
		s2 = ((head_prod - bv) + (head_sumB - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sumB + tail_prod;
		bv = t1 - tail_sumB;
		t2 = ((tail_prod - bv) + (tail_sumB - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sumB = t1 + t2;
		tail_sumB = t2 - (head_sumB - t1);
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumA * split;
	      a11 = con - head_sumA;
	      a11 = con - a11;
	      a21 = head_sumA - a11;
	      con = alpha_i * split;
	      b1 = con - alpha_i;
	      b1 = con - b1;
	      b2 = alpha_i - b1;

	      c11 = head_sumA * alpha_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumA * alpha_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    {
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sumB * split;
	      a11 = con - head_sumB;
	      a11 = con - a11;
	      a21 = head_sumB - a11;
	      con = beta_i * split;
	      b1 = con - beta_i;
	      b1 = con - b1;
	      b2 = beta_i - b1;

	      c11 = head_sumB * beta_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sumB * beta_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp2 = t1 + t2;
	      tail_tmp2 = t2 - (head_tmp2 - t1);
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
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[yi] = head_tmp1;
	    ai += incai;
	    bi += incbi;
	  }
	}
      }
      FPU_FIX_STOP;
    }
    break;

  default:
    {
      BLAS_error(routine_name, -14, prec, 0);
    }
    break;
  }

}				/* end BLAS_dge_sum_mv */
