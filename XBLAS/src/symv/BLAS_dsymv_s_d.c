#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dsymv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta,
		    double *y, int incy)

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
 * a       (input) float*
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
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_dsymv_s_d";

  /* Integer Index Variables */
  int i, k;

  int xi, yi;
  int aik, astarti, x_starti, y_starti;

  int incai;
  int incaik, incaik2;

  int n_i;

  /* Input Matrices */
  const float *a_i = a;
  const double *x_i = x;

  /* Output Vector */
  double *y_i = y;

  /* Input Scalars */
  double alpha_i = alpha;
  double beta_i = beta;

  /* Temporary Floating-Point Variables */
  float a_elem;
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



}				/* end BLAS_dsymv_s_d */
