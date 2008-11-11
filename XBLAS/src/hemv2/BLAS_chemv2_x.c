#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_chemv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, const void *y, int incy,
		   enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a complex Hermitian matrix.
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
 * x_head  (input) void*
 *         Vector x_head
 *
 * x_tail  (input) void*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) const void*
 * 
 * y       (input) void*
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
  const char routine_name[] = "BLAS_chemv2_x";
  switch (prec) {

  case blas_prec_single:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const float *a_i = (float *) a;
      const float *x_head_i = (float *) x_head;
      const float *x_tail_i = (float *) x_tail;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float a_elem[2];
      float x_elem[2];
      float y_elem[2];
      float diag_elem;
      float prod1[2];
      float prod2[2];
      float sum1[2];
      float sum2[2];
      float tmp1[2];
      float tmp2[2];
      float tmp3[2];



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
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

      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod1[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod2[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    prod1[0] = x_elem[0] * diag_elem;
	    prod1[1] = x_elem[1] * diag_elem;
	  }
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    prod2[0] = x_elem[0] * diag_elem;
	    prod2[1] = x_elem[1] * diag_elem;
	  }
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod1[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod2[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] = sum1[0] * alpha_i[0] - sum1[1] * alpha_i[1];
	    tmp1[1] = sum1[0] * alpha_i[1] + sum1[1] * alpha_i[0];
	  }

	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	    tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	  }

	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod1[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod2[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    prod1[0] = x_elem[0] * diag_elem;
	    prod1[1] = x_elem[1] * diag_elem;
	  }
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    prod2[0] = x_elem[0] * diag_elem;
	    prod2[1] = x_elem[1] * diag_elem;
	  }
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod1[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] = a_elem[0] * x_elem[0] - a_elem[1] * x_elem[1];
	      prod2[1] = a_elem[0] * x_elem[1] + a_elem[1] * x_elem[0];
	    }

	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] = sum1[0] * alpha_i[0] - sum1[1] * alpha_i[1];
	    tmp1[1] = sum1[0] * alpha_i[1] + sum1[1] * alpha_i[0];
	  }

	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	    tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	  }

	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
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

      const float *a_i = (float *) a;
      const float *x_head_i = (float *) x_head;
      const float *x_tail_i = (float *) x_tail;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float a_elem[2];
      float x_elem[2];
      float y_elem[2];
      float diag_elem;
      double prod1[2];
      double prod2[2];
      double sum1[2];
      double sum2[2];
      double tmp1[2];
      double tmp2[2];
      double tmp3[2];



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
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

      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod1[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod2[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    prod1[0] = (double) x_elem[0] * diag_elem;
	    prod1[1] = (double) x_elem[1] * diag_elem;
	  }
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    prod2[0] = (double) x_elem[0] * diag_elem;
	    prod2[1] = (double) x_elem[1] * diag_elem;
	  }
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod1[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod2[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  sum1[0] = sum1[1] = 0.0;
	  sum2[0] = sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod1[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod2[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    prod1[0] = (double) x_elem[0] * diag_elem;
	    prod1[1] = (double) x_elem[1] * diag_elem;
	  }
	  sum1[0] = sum1[0] + prod1[0];
	  sum1[1] = sum1[1] + prod1[1];
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    prod2[0] = (double) x_elem[0] * diag_elem;
	    prod2[1] = (double) x_elem[1] * diag_elem;
	  }
	  sum2[0] = sum2[0] + prod2[0];
	  sum2[1] = sum2[1] + prod2[1];
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      prod1[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod1[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum1[0] = sum1[0] + prod1[0];
	    sum1[1] = sum1[1] + prod1[1];
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      prod2[0] =
		(double) a_elem[0] * x_elem[0] -
		(double) a_elem[1] * x_elem[1];
	      prod2[1] =
		(double) a_elem[0] * x_elem[1] +
		(double) a_elem[1] * x_elem[0];
	    }
	    sum2[0] = sum2[0] + prod2[0];
	    sum2[1] = sum2[1] + prod2[1];
	  }
	  sum1[0] = sum1[0] + sum2[0];
	  sum1[1] = sum1[1] + sum2[1];
	  {
	    tmp1[0] =
	      (double) sum1[0] * alpha_i[0] - (double) sum1[1] * alpha_i[1];
	    tmp1[1] =
	      (double) sum1[0] * alpha_i[1] + (double) sum1[1] * alpha_i[0];
	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    tmp2[0] =
	      (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	  }
	  tmp3[0] = tmp1[0] + tmp2[0];
	  tmp3[1] = tmp1[1] + tmp2[1];
	  y_i[yi] = tmp3[0];
	  y_i[yi + 1] = tmp3[1];
	}
      }



      break;
    }

  case blas_prec_extra:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const float *a_i = (float *) a;
      const float *x_head_i = (float *) x_head;
      const float *x_tail_i = (float *) x_tail;
      float *y_i = (float *) y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float a_elem[2];
      float x_elem[2];
      float y_elem[2];
      float diag_elem;
      double head_prod1[2], tail_prod1[2];
      double head_prod2[2], tail_prod2[2];
      double head_sum1[2], tail_sum1[2];
      double head_sum2[2], tail_sum2[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      double head_tmp3[2], tail_tmp3[2];

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

      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaij *= 2;
      incaij2 *= 2;
      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);

      FPU_FIX_START;

      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      if (uplo == blas_lower) {
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  head_sum1[0] = head_sum1[1] = tail_sum1[0] = tail_sum1[1] = 0.0;
	  head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] = 0.0;
	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod1[0] = head_e1;
	      tail_prod1[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod1[1] = head_e1;
	      tail_prod1[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
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
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
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
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod2[0] = head_e1;
	      tail_prod2[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod2[1] = head_e1;
	      tail_prod2[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
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
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
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
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    head_prod1[0] = (double) x_elem[0] * diag_elem;
	    tail_prod1[0] = 0.0;
	    head_prod1[1] = (double) x_elem[1] * diag_elem;
	    tail_prod1[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_prod1[0];
	    tail_b = tail_prod1[0];
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
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_prod1[1];
	    tail_b = tail_prod1[1];
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
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    head_prod2[0] = (double) x_elem[0] * diag_elem;
	    tail_prod2[0] = 0.0;
	    head_prod2[1] = (double) x_elem[1] * diag_elem;
	    tail_prod2[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum2[0];
	    tail_a = tail_sum2[0];
	    head_b = head_prod2[0];
	    tail_b = tail_prod2[0];
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
	    head_sum2[0] = head_t;
	    tail_sum2[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum2[1];
	    tail_a = tail_sum2[1];
	    head_b = head_prod2[1];
	    tail_b = tail_prod2[1];
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
	    head_sum2[1] = head_t;
	    tail_sum2[1] = tail_t;
	  }
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod1[0] = head_e1;
	      tail_prod1[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod1[1] = head_e1;
	      tail_prod1[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
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
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
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
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod2[0] = head_e1;
	      tail_prod2[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod2[1] = head_e1;
	      tail_prod2[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
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
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
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
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_sum2[0];
	    tail_b = tail_sum2[0];
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
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_sum2[1];
	    tail_b = tail_sum2[1];
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
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  {
	    double cd[2];
	    cd[0] = (double) alpha_i[0];
	    cd[1] = (double) alpha_i[1];
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sum1[0];
	      tail_a0 = tail_sum1[0];
	      head_a1 = head_sum1[1];
	      tail_a1 = tail_sum1[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = cd[0] * split;
		b1 = con - cd[0];
		b1 = con - b1;
		b2 = cd[0] - b1;

		c11 = head_a0 * cd[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * cd[0];
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
		con = cd[1] * split;
		b1 = con - cd[1];
		b1 = con - b1;
		b2 = cd[1] - b1;

		c11 = head_a1 * cd[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * cd[1];
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
		con = cd[0] * split;
		b1 = con - cd[0];
		b1 = con - b1;
		b2 = cd[0] - b1;

		c11 = head_a1 * cd[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * cd[0];
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
		con = cd[1] * split;
		b1 = con - cd[1];
		b1 = con - b1;
		b2 = cd[1] - b1;

		c11 = head_a0 * cd[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * cd[1];
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

	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    double head_e1, tail_e1;
	    double d1;
	    double d2;
	    /* Real part */
	    d1 = (double) y_elem[0] * beta_i[0];
	    d2 = (double) -y_elem[1] * beta_i[1];
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
	    head_tmp2[0] = head_e1;
	    tail_tmp2[0] = tail_e1;
	    /* imaginary part */
	    d1 = (double) y_elem[0] * beta_i[1];
	    d2 = (double) y_elem[1] * beta_i[0];
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
	    head_tmp2[1] = head_e1;
	    tail_tmp2[1] = tail_e1;
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
	    head_tmp3[0] = head_t;
	    tail_tmp3[0] = tail_t;
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
	    head_tmp3[1] = head_t;
	    tail_tmp3[1] = tail_t;
	  }
	  y_i[yi] = head_tmp3[0];
	  y_i[yi + 1] = head_tmp3[1];
	}
      } else {
	/* uplo == blas_upper */
	for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	  head_sum1[0] = head_sum1[1] = tail_sum1[0] = tail_sum1[1] = 0.0;
	  head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] = 0.0;

	  for (j = 0, aij = ai, xi = xi0; j < i;
	       j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod1[0] = head_e1;
	      tail_prod1[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod1[1] = head_e1;
	      tail_prod1[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
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
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
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
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod2[0] = head_e1;
	      tail_prod2[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod2[1] = head_e1;
	      tail_prod2[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
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
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
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
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }

	  diag_elem = a_i[aij];
	  x_elem[0] = x_head_i[xi];
	  x_elem[1] = x_head_i[xi + 1];
	  {
	    head_prod1[0] = (double) x_elem[0] * diag_elem;
	    tail_prod1[0] = 0.0;
	    head_prod1[1] = (double) x_elem[1] * diag_elem;
	    tail_prod1[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_prod1[0];
	    tail_b = tail_prod1[0];
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
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_prod1[1];
	    tail_b = tail_prod1[1];
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
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  x_elem[0] = x_tail_i[xi];
	  x_elem[1] = x_tail_i[xi + 1];
	  {
	    head_prod2[0] = (double) x_elem[0] * diag_elem;
	    tail_prod2[0] = 0.0;
	    head_prod2[1] = (double) x_elem[1] * diag_elem;
	    tail_prod2[1] = 0.0;
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum2[0];
	    tail_a = tail_sum2[0];
	    head_b = head_prod2[0];
	    tail_b = tail_prod2[0];
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
	    head_sum2[0] = head_t;
	    tail_sum2[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum2[1];
	    tail_a = tail_sum2[1];
	    head_b = head_prod2[1];
	    tail_b = tail_prod2[1];
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
	    head_sum2[1] = head_t;
	    tail_sum2[1] = tail_t;
	  }
	  j++;
	  aij += incaij2;
	  xi += incx;

	  for (; j < n; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    x_elem[0] = x_head_i[xi];
	    x_elem[1] = x_head_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod1[0] = head_e1;
	      tail_prod1[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod1[1] = head_e1;
	      tail_prod1[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum1[0];
	      tail_a = tail_sum1[0];
	      head_b = head_prod1[0];
	      tail_b = tail_prod1[0];
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
	      head_sum1[0] = head_t;
	      tail_sum1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum1[1];
	      tail_a = tail_sum1[1];
	      head_b = head_prod1[1];
	      tail_b = tail_prod1[1];
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
	      head_sum1[1] = head_t;
	      tail_sum1[1] = tail_t;
	    }
	    x_elem[0] = x_tail_i[xi];
	    x_elem[1] = x_tail_i[xi + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) a_elem[0] * x_elem[0];
	      d2 = (double) -a_elem[1] * x_elem[1];
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
	      head_prod2[0] = head_e1;
	      tail_prod2[0] = tail_e1;
	      /* imaginary part */
	      d1 = (double) a_elem[0] * x_elem[1];
	      d2 = (double) a_elem[1] * x_elem[0];
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
	      head_prod2[1] = head_e1;
	      tail_prod2[1] = tail_e1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_sum2[0];
	      tail_a = tail_sum2[0];
	      head_b = head_prod2[0];
	      tail_b = tail_prod2[0];
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
	      head_sum2[0] = head_t;
	      tail_sum2[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_sum2[1];
	      tail_a = tail_sum2[1];
	      head_b = head_prod2[1];
	      tail_b = tail_prod2[1];
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
	      head_sum2[1] = head_t;
	      tail_sum2[1] = tail_t;
	    }
	  }
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum1[0];
	    tail_a = tail_sum1[0];
	    head_b = head_sum2[0];
	    tail_b = tail_sum2[0];
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
	    head_sum1[0] = head_t;
	    tail_sum1[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum1[1];
	    tail_a = tail_sum1[1];
	    head_b = head_sum2[1];
	    tail_b = tail_sum2[1];
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
	    head_sum1[1] = head_t;
	    tail_sum1[1] = tail_t;
	  }
	  {
	    double cd[2];
	    cd[0] = (double) alpha_i[0];
	    cd[1] = (double) alpha_i[1];
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sum1[0];
	      tail_a0 = tail_sum1[0];
	      head_a1 = head_sum1[1];
	      tail_a1 = tail_sum1[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = cd[0] * split;
		b1 = con - cd[0];
		b1 = con - b1;
		b2 = cd[0] - b1;

		c11 = head_a0 * cd[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * cd[0];
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
		con = cd[1] * split;
		b1 = con - cd[1];
		b1 = con - b1;
		b2 = cd[1] - b1;

		c11 = head_a1 * cd[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * cd[1];
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
		con = cd[0] * split;
		b1 = con - cd[0];
		b1 = con - b1;
		b2 = cd[0] - b1;

		c11 = head_a1 * cd[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * cd[0];
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
		con = cd[1] * split;
		b1 = con - cd[1];
		b1 = con - b1;
		b2 = cd[1] - b1;

		c11 = head_a0 * cd[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * cd[1];
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

	  }
	  y_elem[0] = y_i[yi];
	  y_elem[1] = y_i[yi + 1];
	  {
	    double head_e1, tail_e1;
	    double d1;
	    double d2;
	    /* Real part */
	    d1 = (double) y_elem[0] * beta_i[0];
	    d2 = (double) -y_elem[1] * beta_i[1];
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
	    head_tmp2[0] = head_e1;
	    tail_tmp2[0] = tail_e1;
	    /* imaginary part */
	    d1 = (double) y_elem[0] * beta_i[1];
	    d2 = (double) y_elem[1] * beta_i[0];
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
	    head_tmp2[1] = head_e1;
	    tail_tmp2[1] = tail_e1;
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
	    head_tmp3[0] = head_t;
	    tail_tmp3[0] = tail_t;
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
	    head_tmp3[1] = head_t;
	    tail_tmp3[1] = tail_t;
	  }
	  y_i[yi] = head_tmp3[0];
	  y_i[yi + 1] = head_tmp3[1];
	}
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_chemv2_x */
