#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zhemv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a Hermitian matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input hermitian matrix A.
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
 * x       (input) double*
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
  static const char routine_name[] = "BLAS_zhemv_z_d";

  /* Integer Index Variables */
  int i, k;

  int xi, yi;
  int aik, astarti, x_starti, y_starti;

  int incai;
  int incaik, incaik2;

  int n_i;

  /* Input Matrices */
  const double *a_i = (double *) a;
  const double *x_i = x;

  /* Output Vector */
  double *y_i = (double *) y;

  /* Input Scalars */
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary Floating-Point Variables */
  double a_elem[2];
  double x_elem;
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
  if (order != blas_colmajor && order != blas_rowmajor) {
    BLAS_error(routine_name, -1, order, NULL);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, NULL);
  }
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
  incai *= 2;
  incaik *= 2;
  incaik2 *= 2;
  if (incx < 0) {
    x_starti = (-n_i + 1) * incx;
  } else {
    x_starti = 0;
  }
  if (incy < 0) {
    y_starti = (-n_i + 1) * incy;
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
  } else {
    /*  determine whether we conjugate in first loop or second loop */

    if (uplo == blas_lower) {
      /*  conjugate second */

      /* Case alpha == 1. */
      if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incaik2) {
	    sum[0] = sum[1] = 0.0;
	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aik];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    k++;
	    aik += incaik2;
	    xi += incx;
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
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
	       i < n_i; i++, yi += incy, astarti += incaik2) {
	    sum[0] = sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aik];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    k++;
	    aik += incaik2;
	    xi += incx;
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
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
	     i < n_i; i++, yi += incy, astarti += incaik2) {
	  sum[0] = sum[1] = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];

	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  a_elem[0] = a_i[aik];
	  x_elem = x_i[xi];
	  prod[0] = x_elem * a_elem[0];
	  prod[1] = 0.0;
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	  k++;
	  aik += incaik2;
	  xi += incx;
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
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
    } else {
      /*  conjugate first loop */

      /* Case alpha == 1. */
      if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti, astarti = 0;
	       i < n_i; i++, yi += incy, astarti += incaik2) {
	    sum[0] = sum[1] = 0.0;
	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aik];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    k++;
	    aik += incaik2;
	    xi += incx;
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
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
	       i < n_i; i++, yi += incy, astarti += incaik2) {
	    sum[0] = sum[1] = 0.0;

	    for (k = 0, aik = astarti, xi = x_starti;
		 k < i; k++, aik += incaik, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aik];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    k++;
	    aik += incaik2;
	    xi += incx;
	    for (; k < n_i; k++, aik += incaik2, xi += incx) {
	      a_elem[0] = a_i[aik];
	      a_elem[1] = a_i[aik + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
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
	     i < n_i; i++, yi += incy, astarti += incaik2) {
	  sum[0] = sum[1] = 0.0;

	  for (k = 0, aik = astarti, xi = x_starti;
	       k < i; k++, aik += incaik, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  a_elem[0] = a_i[aik];
	  x_elem = x_i[xi];
	  prod[0] = x_elem * a_elem[0];
	  prod[1] = 0.0;
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	  k++;
	  aik += incaik2;
	  xi += incx;
	  for (; k < n_i; k++, aik += incaik2, xi += incx) {
	    a_elem[0] = a_i[aik];
	    a_elem[1] = a_i[aik + 1];

	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
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
    }
  }


}				/* end BLAS_zhemv_z_d */
