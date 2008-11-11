#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_zhemv2_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, const void *y, int incy)

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
 */
{
  /* Routine name */
  const char routine_name[] = "BLAS_zhemv2_c_c";

  int i, j;
  int xi, yi, xi0, yi0;
  int aij, ai;
  int incai;
  int incaij, incaij2;

  const float *a_i = (float *) a;
  const float *x_head_i = (float *) x_head;
  const float *x_tail_i = (float *) x_tail;
  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float a_elem[2];
  float x_elem[2];
  double y_elem[2];
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
      for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	a_elem[0] = a_i[aij];
	a_elem[1] = a_i[aij + 1];
	x_elem[0] = x_head_i[xi];
	x_elem[1] = x_head_i[xi + 1];
	{
	  prod1[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod1[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	}
	sum1[0] = sum1[0] + prod1[0];
	sum1[1] = sum1[1] + prod1[1];
	x_elem[0] = x_tail_i[xi];
	x_elem[1] = x_tail_i[xi + 1];
	{
	  prod2[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod2[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
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
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod1[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	}
	sum1[0] = sum1[0] + prod1[0];
	sum1[1] = sum1[1] + prod1[1];
	x_elem[0] = x_tail_i[xi];
	x_elem[1] = x_tail_i[xi + 1];
	{
	  prod2[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod2[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
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

      for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	a_elem[0] = a_i[aij];
	a_elem[1] = a_i[aij + 1];
	a_elem[1] = -a_elem[1];
	x_elem[0] = x_head_i[xi];
	x_elem[1] = x_head_i[xi + 1];
	{
	  prod1[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod1[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	}
	sum1[0] = sum1[0] + prod1[0];
	sum1[1] = sum1[1] + prod1[1];
	x_elem[0] = x_tail_i[xi];
	x_elem[1] = x_tail_i[xi + 1];
	{
	  prod2[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod2[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
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
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod1[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	}
	sum1[0] = sum1[0] + prod1[0];
	sum1[1] = sum1[1] + prod1[1];
	x_elem[0] = x_tail_i[xi];
	x_elem[1] = x_tail_i[xi + 1];
	{
	  prod2[0] =
	    (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	  prod2[1] =
	    (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
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



}				/* end BLAS_zhemv2_c_c */
