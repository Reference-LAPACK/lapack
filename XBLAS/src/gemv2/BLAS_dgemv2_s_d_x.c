#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dgemv2_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, double alpha, const float *a, int lda,
		       const double *head_x, const double *tail_x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * op(A) * head_x + alpha * op(A) * tail_x + beta * y,
 * where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of A: no trans, trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of A
 *
 * n            (input) int
 *              Dimension of A and the length of vector x and z
 *
 * alpha        (input) double
 *              
 * A            (input) const float*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * head_x
 * tail_x       (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) double
 *
 * y            (input) const double*
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
  static const char routine_name[] = "BLAS_dgemv2_s_d_x";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const float *a_i = a;
      const double *head_x_i = head_x;
      const double *tail_x_i = tail_x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float a_elem;
      double x_elem;
      double y_elem;
      double prod;
      double sum;
      double sum2;
      double tmp1;
      double tmp2;


      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -10, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }

      if (lda < leny)
	BLAS_error(routine_name, -7, lda, NULL);








      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem = y_i[iy];
	    tmp1 = y_elem * beta_i;
	    y_i[iy] = tmp1;
	    iy += incy;
	  }
	}
      } else {			/* alpha != 0 */

	/* if beta = 0, we can save m multiplies:
	   y = alpha*A*head_x + alpha*A*tail_x  */
	if (beta_i == 0.0) {
	  if (alpha_i == 1.0) {
	    /* save m more multiplies if alpha = 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      sum2 = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		x_elem = tail_x_i[jx];
		prod = a_elem * x_elem;
		sum2 = sum2 + prod;
		aij += incaij;
		jx += incx;
	      }
	      sum = sum + sum2;
	      y_i[iy] = sum;
	      ai += incai;
	      iy += incy;
	    }			/* end for */
	  } else {		/* alpha != 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      sum2 = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		x_elem = tail_x_i[jx];
		prod = a_elem * x_elem;
		sum2 = sum2 + prod;
		aij += incaij;
		jx += incx;
	      }
	      tmp1 = sum * alpha_i;
	      tmp2 = sum2 * alpha_i;
	      tmp1 = tmp1 + tmp2;
	      y_i[iy] = tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {		/* beta != 0 */
	  if (alpha_i == 1.0) {
	    /* save m multiplies if alpha = 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;;
	      sum2 = 0.0;;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		x_elem = tail_x_i[jx];
		prod = a_elem * x_elem;
		sum2 = sum2 + prod;
		aij += incaij;
		jx += incx;
	      }
	      sum = sum + sum2;
	      y_elem = y_i[iy];
	      tmp1 = y_elem * beta_i;
	      tmp2 = sum + tmp1;
	      y_i[iy] = tmp2;
	      ai += incai;
	      iy += incy;
	    }
	  } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;;
	      sum2 = 0.0;;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		x_elem = tail_x_i[jx];
		prod = a_elem * x_elem;
		sum2 = sum2 + prod;
		aij += incaij;
		jx += incx;
	      }
	      tmp1 = sum * alpha_i;
	      tmp2 = sum2 * alpha_i;
	      tmp1 = tmp1 + tmp2;
	      y_elem = y_i[iy];
	      tmp2 = y_elem * beta_i;
	      tmp1 = tmp1 + tmp2;
	      y_i[iy] = tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	}

      }



      break;
    }
  case blas_prec_extra:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const float *a_i = a;
      const double *head_x_i = head_x;
      const double *tail_x_i = tail_x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float a_elem;
      double x_elem;
      double y_elem;
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_sum2, tail_sum2;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;
      FPU_FIX_DECL;

      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -10, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }

      if (lda < leny)
	BLAS_error(routine_name, -7, lda, NULL);

      FPU_FIX_START;






      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem = y_i[iy];
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
	    y_i[iy] = head_tmp1;
	    iy += incy;
	  }
	}
      } else {			/* alpha != 0 */

	/* if beta = 0, we can save m multiplies:
	   y = alpha*A*head_x + alpha*A*tail_x  */
	if (beta_i == 0.0) {
	  if (alpha_i == 1.0) {
	    /* save m more multiplies if alpha = 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      head_sum2 = tail_sum2 = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
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
		x_elem = tail_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum2 + head_prod;
		  bv = s1 - head_sum2;
		  s2 = ((head_prod - bv) + (head_sum2 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum2 + tail_prod;
		  bv = t1 - tail_sum2;
		  t2 = ((tail_prod - bv) + (tail_sum2 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum2 = t1 + t2;
		  tail_sum2 = t2 - (head_sum2 - t1);
		}
		aij += incaij;
		jx += incx;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_sum2;
		bv = s1 - head_sum;
		s2 = ((head_sum2 - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_sum2;
		bv = t1 - tail_sum;
		t2 = ((tail_sum2 - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	      y_i[iy] = head_sum;
	      ai += incai;
	      iy += incy;
	    }			/* end for */
	  } else {		/* alpha != 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      head_sum2 = tail_sum2 = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
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
		x_elem = tail_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum2 + head_prod;
		  bv = s1 - head_sum2;
		  s2 = ((head_prod - bv) + (head_sum2 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum2 + tail_prod;
		  bv = t1 - tail_sum2;
		  t2 = ((tail_prod - bv) + (tail_sum2 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum2 = t1 + t2;
		  tail_sum2 = t2 - (head_sum2 - t1);
		}
		aij += incaij;
		jx += incx;
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
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum2 * split;
		a11 = con - head_sum2;
		a11 = con - a11;
		a21 = head_sum2 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_sum2 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum2 * alpha_i;
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
	      y_i[iy] = head_tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {		/* beta != 0 */
	  if (alpha_i == 1.0) {
	    /* save m multiplies if alpha = 1 */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;;
	      head_sum2 = tail_sum2 = 0.0;;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
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
		x_elem = tail_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum2 + head_prod;
		  bv = s1 - head_sum2;
		  s2 = ((head_prod - bv) + (head_sum2 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum2 + tail_prod;
		  bv = t1 - tail_sum2;
		  t2 = ((tail_prod - bv) + (tail_sum2 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum2 = t1 + t2;
		  tail_sum2 = t2 - (head_sum2 - t1);
		}
		aij += incaij;
		jx += incx;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_sum2;
		bv = s1 - head_sum;
		s2 = ((head_sum2 - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_sum2;
		bv = t1 - tail_sum;
		t2 = ((tail_sum2 - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	      y_elem = y_i[iy];
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
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_tmp1;
		bv = s1 - head_sum;
		s2 = ((head_tmp1 - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_tmp1;
		bv = t1 - tail_sum;
		t2 = ((tail_tmp1 - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_tmp2 = t1 + t2;
		tail_tmp2 = t2 - (head_tmp2 - t1);
	      }
	      y_i[iy] = head_tmp2;
	      ai += incai;
	      iy += incy;
	    }
	  } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;;
	      head_sum2 = tail_sum2 = 0.0;;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = head_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
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
		x_elem = tail_x_i[jx];
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

		    head_prod = dt * x_elem;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum2 + head_prod;
		  bv = s1 - head_sum2;
		  s2 = ((head_prod - bv) + (head_sum2 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum2 + tail_prod;
		  bv = t1 - tail_sum2;
		  t2 = ((tail_prod - bv) + (tail_sum2 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum2 = t1 + t2;
		  tail_sum2 = t2 - (head_sum2 - t1);
		}
		aij += incaij;
		jx += incx;
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
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum2 * split;
		a11 = con - head_sum2;
		a11 = con - a11;
		a21 = head_sum2 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_sum2 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum2 * alpha_i;
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
	      y_elem = y_i[iy];
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
	      y_i[iy] = head_tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	}

      }

      FPU_FIX_STOP;
    }
    break;
  }
}
