#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, const void *alpha,
		  const void *a, int lda, const void *x, int incx,
		  const void *beta, void *y, int incy,
		  enum blas_prec_type prec)

/*           
 * Purpose
 * =======
 *
 *  gbmv computes y = alpha * A * x + beta * y, where 
 *
 *  A is a m x n banded matrix
 *  x is a n x 1 vector
 *  y is a m x 1 vector
 *  alpha and beta are scalars 
 *
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
 * kl           (input) int 
 *              Number of lower diagnols of AB
 *
 * ku           (input) int
 *              Number of upper diagnols of AB
 *
 * alpha        (input) const void*
 *              
 * AB           (input) void*
 *
 * lda          (input) int 
 *              Leading dimension of AB
 *              lda >= ku + kl + 1
 *
 * x            (input) void*
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
 *
 * LOCAL VARIABLES 
 * ===============
 * 
 *  As an example, these variables are described on the mxn, column 
 *  major, banded matrix described in section 2.2.3 of the specification  
 *
 *  astart      indexes first element in A where computation begins
 *
 *  incai1      indexes first element in row where row is less than lbound
 * 
 *  incai2      indexes first element in row where row exceeds lbound
 *   
 *  lbound      denotes the number of rows before  first element shifts 
 *
 *  rbound      denotes the columns where there is blank space
 *   
 *  ra          index of the rightmost element for a given row
 *  
 *  la          index of leftmost  elements for a given row
 *
 *  ra - la     width of a row
 *
 *                        rbound 
 *            la   ra    ____|_____ 
 *             |    |   |          |
 *         |  a00  a01   *    *   *
 * lbound -|  a10  a11  a12   *   *
 *         |  a20  a21  a22  a23  *
 *             *   a31  a32  a33 a34
 *             *    *   a42  a43 a44
 *
 *  Varations on order and transpose have been implemented by modifying these
 *  local variables. 
 *
 */
{
  static const char routine_name[] = "BLAS_zgbmv_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int ky, iy, kx, jx, j, i, rbound, lbound, ra, la, lenx, leny;
      int incaij, aij, incai1, incai2, astart, ai;
      double *y_i = (double *) y;
      const double *a_i = (double *) a;
      const double *x_i = (double *) x;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double tmp1[2];
      double tmp2[2];
      double result[2];
      double sum[2];
      double prod[2];
      double a_elem[2];
      double x_elem[2];
      double y_elem[2];


      if (order != blas_colmajor && order != blas_rowmajor)
	BLAS_error(routine_name, -1, order, NULL);
      if (trans != blas_no_trans &&
	  trans != blas_trans && trans != blas_conj_trans) {
	BLAS_error(routine_name, -2, trans, NULL);
      }
      if (m < 0)
	BLAS_error(routine_name, -3, m, NULL);
      if (n < 0)
	BLAS_error(routine_name, -4, n, NULL);
      if (kl < 0 || kl >= m)
	BLAS_error(routine_name, -5, kl, NULL);
      if (ku < 0 || ku >= n)
	BLAS_error(routine_name, -6, ku, NULL);
      if (lda < kl + ku + 1)
	BLAS_error(routine_name, -9, lda, NULL);
      if (incx == 0)
	BLAS_error(routine_name, -11, incx, NULL);
      if (incy == 0)
	BLAS_error(routine_name, -14, incy, NULL);

      if ((m == 0) || (n == 0) ||
	  (((alpha_i[0] == 0.0 && alpha_i[1] == 0.0)
	    && ((beta_i[0] == 1.0 && beta_i[1] == 0.0)))))
	return;

      if (trans == blas_no_trans) {
	lenx = n;
	leny = m;
      } else {			/* change back */
	lenx = m;
	leny = n;
      }

      if (incx < 0) {
	kx = -(lenx - 1) * incx;
      } else {
	kx = 0;
      }

      if (incy < 0) {
	ky = -(leny - 1) * incy;
      } else {
	ky = 0;
      }



      /* if alpha = 0, return y = y*beta */
      if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	astart = ku;
	incai1 = 1;
	incai2 = lda;
	incaij = lda - 1;
	lbound = kl;
	rbound = n - ku - 1;
	ra = ku;
      } else if ((order == blas_colmajor) && (trans != blas_no_trans)) {
	astart = ku;
	incai1 = lda - 1;
	incai2 = lda;
	incaij = 1;
	lbound = ku;
	rbound = m - kl - 1;
	ra = kl;
      } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	astart = kl;
	incai1 = lda - 1;
	incai2 = lda;
	incaij = 1;
	lbound = kl;
	rbound = n - ku - 1;
	ra = ku;
      } else {			/* rowmajor and blas_trans */
	astart = kl;
	incai1 = 1;
	incai2 = lda;
	incaij = lda - 1;
	lbound = ku;
	rbound = m - kl - 1;
	ra = kl;
      }
      incx *= 2;
      incy *= 2;
      incaij *= 2;
      incai1 *= 2;
      incai2 *= 2;
      astart *= 2;
      ky *= 2;
      kx *= 2;

      la = 0;
      ai = astart;
      iy = ky;
      for (i = 0; i < leny; i++) {
	sum[0] = sum[1] = 0.0;
	aij = ai;
	jx = kx;
	if (trans != blas_conj_trans) {
	  for (j = ra - la; j >= 0; j--) {
	    x_elem[0] = x_i[jx];
	    x_elem[1] = x_i[jx + 1];
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    {
	      prod[0] =
		(double) x_elem[0] * a_elem[0] -
		(double) x_elem[1] * a_elem[1];
	      prod[1] =
		(double) x_elem[0] * a_elem[1] +
		(double) x_elem[1] * a_elem[0];
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij;
	    jx += incx;
	  }

	} else {
	  for (j = ra - la; j >= 0; j--) {
	    x_elem[0] = x_i[jx];
	    x_elem[1] = x_i[jx + 1];
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    {
	      prod[0] =
		(double) x_elem[0] * a_elem[0] -
		(double) x_elem[1] * a_elem[1];
	      prod[1] =
		(double) x_elem[0] * a_elem[1] +
		(double) x_elem[1] * a_elem[0];
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij;
	    jx += incx;
	  }
	}

	{
	  tmp1[0] =
	    (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	  tmp1[1] =
	    (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	}
	y_elem[0] = y_i[iy];
	y_elem[1] = y_i[iy + 1];
	{
	  tmp2[0] =
	    (double) beta_i[0] * y_elem[0] - (double) beta_i[1] * y_elem[1];
	  tmp2[1] =
	    (double) beta_i[0] * y_elem[1] + (double) beta_i[1] * y_elem[0];
	}
	result[0] = tmp1[0] + tmp2[0];
	result[1] = tmp1[1] + tmp2[1];
	y_i[iy] = result[0];
	y_i[iy + 1] = result[1];
	iy += incy;
	if (i >= lbound) {
	  kx += incx;
	  ai += incai2;
	  la++;
	} else {
	  ai += incai1;
	}
	if (i < rbound) {
	  ra++;
	}
      }


    }
    break;
  case blas_prec_extra:
    {
      int ky, iy, kx, jx, j, i, rbound, lbound, ra, la, lenx, leny;
      int incaij, aij, incai1, incai2, astart, ai;
      double *y_i = (double *) y;
      const double *a_i = (double *) a;
      const double *x_i = (double *) x;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      double result[2];
      double head_sum[2], tail_sum[2];
      double head_prod[2], tail_prod[2];
      double a_elem[2];
      double x_elem[2];
      double y_elem[2];
      FPU_FIX_DECL;

      if (order != blas_colmajor && order != blas_rowmajor)
	BLAS_error(routine_name, -1, order, NULL);
      if (trans != blas_no_trans &&
	  trans != blas_trans && trans != blas_conj_trans) {
	BLAS_error(routine_name, -2, trans, NULL);
      }
      if (m < 0)
	BLAS_error(routine_name, -3, m, NULL);
      if (n < 0)
	BLAS_error(routine_name, -4, n, NULL);
      if (kl < 0 || kl >= m)
	BLAS_error(routine_name, -5, kl, NULL);
      if (ku < 0 || ku >= n)
	BLAS_error(routine_name, -6, ku, NULL);
      if (lda < kl + ku + 1)
	BLAS_error(routine_name, -9, lda, NULL);
      if (incx == 0)
	BLAS_error(routine_name, -11, incx, NULL);
      if (incy == 0)
	BLAS_error(routine_name, -14, incy, NULL);

      if ((m == 0) || (n == 0) ||
	  (((alpha_i[0] == 0.0 && alpha_i[1] == 0.0)
	    && ((beta_i[0] == 1.0 && beta_i[1] == 0.0)))))
	return;

      if (trans == blas_no_trans) {
	lenx = n;
	leny = m;
      } else {			/* change back */
	lenx = m;
	leny = n;
      }

      if (incx < 0) {
	kx = -(lenx - 1) * incx;
      } else {
	kx = 0;
      }

      if (incy < 0) {
	ky = -(leny - 1) * incy;
      } else {
	ky = 0;
      }

      FPU_FIX_START;

      /* if alpha = 0, return y = y*beta */
      if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	astart = ku;
	incai1 = 1;
	incai2 = lda;
	incaij = lda - 1;
	lbound = kl;
	rbound = n - ku - 1;
	ra = ku;
      } else if ((order == blas_colmajor) && (trans != blas_no_trans)) {
	astart = ku;
	incai1 = lda - 1;
	incai2 = lda;
	incaij = 1;
	lbound = ku;
	rbound = m - kl - 1;
	ra = kl;
      } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	astart = kl;
	incai1 = lda - 1;
	incai2 = lda;
	incaij = 1;
	lbound = kl;
	rbound = n - ku - 1;
	ra = ku;
      } else {			/* rowmajor and blas_trans */
	astart = kl;
	incai1 = 1;
	incai2 = lda;
	incaij = lda - 1;
	lbound = ku;
	rbound = m - kl - 1;
	ra = kl;
      }
      incx *= 2;
      incy *= 2;
      incaij *= 2;
      incai1 *= 2;
      incai2 *= 2;
      astart *= 2;
      ky *= 2;
      kx *= 2;

      la = 0;
      ai = astart;
      iy = ky;
      for (i = 0; i < leny; i++) {
	head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
	aij = ai;
	jx = kx;
	if (trans != blas_conj_trans) {
	  for (j = ra - la; j >= 0; j--) {
	    x_elem[0] = x_i[jx];
	    x_elem[1] = x_i[jx + 1];
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
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
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t1 = x_elem[0] * a_elem[0];
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
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t2 = x_elem[1] * a_elem[1];
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
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t1 = x_elem[1] * a_elem[0];
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
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t2 = x_elem[0] * a_elem[1];
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
	    aij += incaij;
	    jx += incx;
	  }

	} else {
	  for (j = ra - la; j >= 0; j--) {
	    x_elem[0] = x_i[jx];
	    x_elem[1] = x_i[jx + 1];
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
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
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t1 = x_elem[0] * a_elem[0];
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
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t2 = x_elem[1] * a_elem[1];
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
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t1 = x_elem[1] * a_elem[0];
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
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t2 = x_elem[0] * a_elem[1];
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
	    aij += incaij;
	    jx += incx;
	  }
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

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_elem[0] * split;
	    b1 = con - y_elem[0];
	    b1 = con - b1;
	    b2 = y_elem[0] - b1;

	    head_t1 = beta_i[0] * y_elem[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_elem[1] * split;
	    b1 = con - y_elem[1];
	    b1 = con - b1;
	    b2 = y_elem[1] - b1;

	    head_t2 = beta_i[1] * y_elem[1];
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

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_elem[0] * split;
	    b1 = con - y_elem[0];
	    b1 = con - b1;
	    b2 = y_elem[0] - b1;

	    head_t1 = beta_i[1] * y_elem[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_elem[1] * split;
	    b1 = con - y_elem[1];
	    b1 = con - b1;
	    b2 = y_elem[1] - b1;

	    head_t2 = beta_i[0] * y_elem[1];
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
	    result[0] = t1 + t2;
	  }
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
	    result[1] = t1 + t2;
	  }
	}
	y_i[iy] = result[0];
	y_i[iy + 1] = result[1];
	iy += incy;
	if (i >= lbound) {
	  kx += incx;
	  ai += incai2;
	  la++;
	} else {
	  ai += incai1;
	}
	if (i < rbound) {
	  ra++;
	}
      }

      FPU_FIX_STOP;
    }
    break;
  }
}				/* end GEMV_NAME(z, z, z, _x) */
