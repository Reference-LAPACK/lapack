#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_csymv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy)

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
 * a       (input) float*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x       (input) float*
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
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_csymv_s_s";

  /* Integer Index Variables */
  int i, k;

  int xi, yi;
  int aik, astarti, x_starti, y_starti;

  int incai;
  int incaik, incaik2;

  int n_i;

  /* Input Matrices */
  const float *a_i = a;
  const float *x_i = x;

  /* Output Vector */
  float *y_i = (float *) y;

  /* Input Scalars */
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary Floating-Point Variables */
  float a_elem;
  float x_elem;
  float y_elem[2];
  float prod;
  float sum;
  float tmp1[2];
  float tmp2[2];



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

  incy *= 2;



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
	tmp1[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	tmp1[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
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
	tmp1[0] = sum;
	tmp1[1] = 0.0;
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
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
	y_elem[0] = y_i[yi];
	y_elem[1] = y_i[yi + 1];
	{
	  tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	  tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	}

	tmp1[0] = sum;
	tmp1[1] = 0.0;
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
      y_elem[0] = y_i[yi];
      y_elem[1] = y_i[yi + 1];
      {
	tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
      }

      {
	tmp1[0] = alpha_i[0] * sum;
	tmp1[1] = alpha_i[1] * sum;
      }
      tmp1[0] = tmp2[0] + tmp1[0];
      tmp1[1] = tmp2[1] + tmp1[1];
      y_i[yi] = tmp1[0];
      y_i[yi + 1] = tmp1[1];
    }
  }



}				/* end BLAS_csymv_s_s */
