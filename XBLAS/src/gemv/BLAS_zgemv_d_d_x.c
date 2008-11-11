#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zgemv_d_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const double *a,
		      int lda, const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * y, where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of AP; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of AB; no trans, 
 *              trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of AB
 *
 * n            (input) int
 *              Dimension of AB and the length of vector x
 *
 * alpha        (input) const void*
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
 * beta         (input) const void*
 *
 * y            (input/output) void*
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
  static const char routine_name[] = "BLAS_zgemv_d_d_x";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const double *a_i = a;
      const double *x_i = x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem;
      double x_elem;
      double y_elem[2];
      double prod;
      double sum;
      double tmp1[2];
      double tmp2[2];


      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -9, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -12, incy, 0);

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
      if ((order == blas_colmajor && lda < m) ||
	  (order == blas_rowmajor && lda < n))
	BLAS_error(routine_name, -7, lda, NULL);




      incy *= 2;



      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    y_i[iy + 1] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i[0] == 0.0 && beta_i[1] == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
	    {
	      tmp1[0] =
		(double) y_elem[0] * beta_i[0] -
		(double) y_elem[1] * beta_i[1];
	      tmp1[1] =
		(double) y_elem[0] * beta_i[1] +
		(double) y_elem[1] * beta_i[0];
	    }
	    y_i[iy] = tmp1[0];
	    y_i[iy + 1] = tmp1[1];
	    iy += incy;
	  }
	}
      } else {

	/* if beta = 0, we can save m multiplies: y = alpha*A*x */
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* save m more multiplies if alpha = 1 */
	  if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		aij += incaij;
		jx += incx;
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      y_i[iy] = tmp1[0];
	      y_i[iy + 1] = tmp1[1];
	      ai += incai;
	      iy += incy;
	    }
	  } else {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		aij += incaij;
		jx += incx;
	      }
	      {
		tmp1[0] = alpha_i[0] * sum;
		tmp1[1] = alpha_i[1] * sum;
	      }
	      y_i[iy] = tmp1[0];
	      y_i[iy + 1] = tmp1[1];
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {
	  /* the most general form, y = alpha*A*x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem = a_i[aij];

	      x_elem = x_i[jx];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      tmp1[0] = alpha_i[0] * sum;
	      tmp1[1] = alpha_i[1] * sum;
	    }
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
	    {
	      tmp2[0] =
		(double) y_elem[0] * beta_i[0] -
		(double) y_elem[1] * beta_i[1];
	      tmp2[1] =
		(double) y_elem[0] * beta_i[1] +
		(double) y_elem[1] * beta_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[iy] = tmp1[0];
	    y_i[iy + 1] = tmp1[1];
	    ai += incai;
	    iy += incy;
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

      const double *a_i = a;
      const double *x_i = x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem;
      double x_elem;
      double y_elem[2];
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      FPU_FIX_DECL;

      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -9, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -12, incy, 0);

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
      if ((order == blas_colmajor && lda < m) ||
	  (order == blas_rowmajor && lda < n))
	BLAS_error(routine_name, -7, lda, NULL);

      FPU_FIX_START;


      incy *= 2;



      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    y_i[iy + 1] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i[0] == 0.0 && beta_i[1] == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
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
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }
	    y_i[iy] = head_tmp1[0];
	    y_i[iy + 1] = head_tmp1[1];
	    iy += incy;
	  }
	}
      } else {

	/* if beta = 0, we can save m multiplies: y = alpha*A*x */
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* save m more multiplies if alpha = 1 */
	  if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
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
		aij += incaij;
		jx += incx;
	      }
	      head_tmp1[0] = head_sum;
	      tail_tmp1[0] = tail_sum;
	      head_tmp1[1] = tail_tmp1[1] = 0.0;
	      y_i[iy] = head_tmp1[0];
	      y_i[iy + 1] = head_tmp1[1];
	      ai += incai;
	      iy += incy;
	    }
	  } else {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
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
		aij += incaij;
		jx += incx;
	      }
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = alpha_i[0] * split;
		  b1 = con - alpha_i[0];
		  b1 = con - b1;
		  b2 = alpha_i[0] - b1;

		  c11 = head_sum * alpha_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * alpha_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_tmp1[0] = head_t;
		tail_tmp1[0] = tail_t;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = alpha_i[1] * split;
		  b1 = con - alpha_i[1];
		  b1 = con - b1;
		  b2 = alpha_i[1] - b1;

		  c11 = head_sum * alpha_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * alpha_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_tmp1[1] = head_t;
		tail_tmp1[1] = tail_t;
	      }
	      y_i[iy] = head_tmp1[0];
	      y_i[iy + 1] = head_tmp1[1];
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {
	  /* the most general form, y = alpha*A*x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    head_sum = tail_sum = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem = a_i[aij];

	      x_elem = x_i[jx];
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
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_sum * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_sum * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
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
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
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
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
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
	    y_i[iy] = head_tmp1[0];
	    y_i[iy + 1] = head_tmp1[1];
	    ai += incai;
	    iy += incy;
	  }
	}

      }

      FPU_FIX_STOP;
    }
    break;
  }
}
