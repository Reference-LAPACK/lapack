#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_chbmv_c_s(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
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
 * where A is a hermitian band matrix.
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
 * k       (input) int
 *         Number of subdiagonals ( = number of superdiagonals)
 *      
 * alpha   (input) const void*
 * 
 * a       (input) void*
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
 *
 *  Notes on storing a hermitian band matrix:
 * 
 *      Integers in the below arrays represent values of
 *              type complex float.
 *
 *    if we have a hermitian matrix:
 *
 *      1d  2   3   0   0
 *      2#  4d  5   6   0
 *      3#  5#  7d  8   9
 *      0   6#  8#  10d 11
 *      0   0   9#  11# 12d
 *
 *     This matrix has n == 5, and k == 2. It can be stored in the
 *      following ways:
 *
 *      Notes for the examples:
 *      Each column below represents a contigous vector.
 *      Columns are strided by lda.
 *      An asterisk (*) represents a position in the 
 *       matrix that is not used.
 *      A pound sign (#) represents the conjugated form is stored
 *      A d following an integer indicates that the imaginary
 *       part of the number is assumed to be zero.
 *      Note that the minimum lda (size of column) is 3 (k+1).
 *       lda may be arbitrarily large; an lda > 3 would mean
 *       there would be unused data at the bottom of the below
 *       columns.        
 *
 *    blas_colmajor and blas_upper:
 *      *   *   3   6   9
 *      *   2   5   8   11 
 *      1d  4d  7d  10d 12d
 *
 *
 *    blas_colmajor and blas_lower
 *      1d   4d   7d   10d  12d
 *      2#   5#   8#   11#  *
 *      3#   6#   9#   *    *
 *
 *
 *    blas_rowmajor and blas_upper 
 *      Columns here also represent contigous arrays.
 *      1d  4d  7d  10d  12d
 *      2   5   8   11   *
 *      3   6   9   *    *
 *
 *
 *    blas_rowmajor and blas_lower
 *      Columns here also represent contigous arrays.
 *      *   *   3#  6#   9#
 *      *   2#  5#  8#   11#
 *      1d  4d  7d  10d  12d
 *
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_chbmv_c_s";

  /* Integer Index Variables */
  int i, j;
  int xi, yi;
  int aij, astarti, x_starti, y_starti;
  int incaij, incaij2;
  int n_i;
  int maxj_first, maxj_second;

  /* Input Matrices */
  const float *a_i = (float *) a;
  const float *x_i = x;

  /* Output Vector */
  float *y_i = (float *) y;

  /* Input Scalars */
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary Floating-Point Variables */
  float a_elem[2];
  float x_elem;
  float y_elem[2];
  float prod[2];
  float sum[2];
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
  if (order != blas_colmajor && order != blas_rowmajor) {
    BLAS_error(routine_name, -1, order, 0);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, 0);
  }
  if (n < 0) {
    BLAS_error(routine_name, -3, n, 0);
  }
  if (k < 0 || k > n) {
    BLAS_error(routine_name, -4, k, 0);
  }
  if ((lda < k + 1) || (lda < 1)) {
    BLAS_error(routine_name, -7, lda, 0);
  }
  if (incx == 0) {
    BLAS_error(routine_name, -9, incx, 0);
  }
  if (incy == 0) {
    BLAS_error(routine_name, -12, incy, 0);
  }

  /* Set Index Parameters */
  n_i = n;

  if (((uplo == blas_upper) && (order == blas_colmajor)) ||
      ((uplo == blas_lower) && (order == blas_rowmajor))) {
    incaij = 1;			/* increment in first loop */
    incaij2 = lda - 1;		/* increment in second loop */
    astarti = k;		/* does not start on zero element */
  } else {
    incaij = lda - 1;
    incaij2 = 1;
    astarti = 0;		/* start on first element of array */
  }
  /* Adjustment to increments (if any) */
  incy *= 2;
  astarti *= 2;
  incaij *= 2;
  incaij2 *= 2;
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
	tmp1[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	tmp1[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
      }

      y_i[yi] = tmp1[0];
      y_i[yi + 1] = tmp1[1];
    }
  } else {
    /*  determine the loop interation counts */
    /* maj_first is number of elements done in first loop 
       (this will increase by one over each column up to a limit) */
    maxj_first = 0;

    /* maxj_second is number of elements done in 
       second loop the first time */
    maxj_second = MIN(k + 1, n_i);

    /*  determine whether we conjugate in first loop or second loop */
    if (uplo == blas_lower) {
      /*  conjugate second loop */

      /* Case alpha == 1. */
      if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	    sum[0] = sum[1] = 0.0;
	    for (j = 0, aij = astarti, xi = x_starti;
		 j < maxj_first; j++, aij += incaij, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aij];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij2;
	    xi += incx;
	    for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
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
	    if (i + 1 >= (n_i - k)) {
	      maxj_second--;
	    }
	    if (i >= k) {
	      astarti += (incaij + incaij2);
	      x_starti += incx;
	    } else {
	      maxj_first++;
	      astarti += incaij2;
	    }
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	    sum[0] = sum[1] = 0.0;

	    for (j = 0, aij = astarti, xi = x_starti;
		 j < maxj_first; j++, aij += incaij, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aij];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij2;
	    xi += incx;
	    for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
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
	      tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	      tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	    }

	    tmp1[0] = sum[0];
	    tmp1[1] = sum[1];
	    tmp1[0] = tmp2[0] + tmp1[0];
	    tmp1[1] = tmp2[1] + tmp1[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    if (i + 1 >= (n_i - k)) {
	      maxj_second--;
	    }
	    if (i >= k) {
	      astarti += (incaij + incaij2);
	      x_starti += incx;
	    } else {
	      maxj_first++;
	      astarti += incaij2;
	    }
	  }
	}
      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum[0] = sum[1] = 0.0;

	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];

	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  a_elem[0] = a_i[aij];
	  x_elem = x_i[xi];
	  prod[0] = x_elem * a_elem[0];
	  prod[1] = 0.0;
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	  aij += incaij2;
	  xi += incx;
	  for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
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
	    tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	    tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	  }

	  {
	    tmp1[0] = sum[0] * alpha_i[0] - sum[1] * alpha_i[1];
	    tmp1[1] = sum[0] * alpha_i[1] + sum[1] * alpha_i[0];
	  }

	  tmp1[0] = tmp2[0] + tmp1[0];
	  tmp1[1] = tmp2[1] + tmp1[1];
	  y_i[yi] = tmp1[0];
	  y_i[yi + 1] = tmp1[1];
	  if (i + 1 >= (n_i - k)) {
	    maxj_second--;
	  }
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      }
    } else {
      /*  conjugate first loop */

      /* Case alpha == 1. */
      if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	  for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	    sum[0] = sum[1] = 0.0;
	    for (j = 0, aij = astarti, xi = x_starti;
		 j < maxj_first; j++, aij += incaij, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aij];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij2;
	    xi += incx;
	    for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

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
	    if (i + 1 >= (n_i - k)) {
	      maxj_second--;
	    }
	    if (i >= k) {
	      astarti += (incaij + incaij2);
	      x_starti += incx;
	    } else {
	      maxj_first++;
	      astarti += incaij2;
	    }
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  y  <--- A * x + beta * y */
	  for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	    sum[0] = sum[1] = 0.0;

	    for (j = 0, aij = astarti, xi = x_starti;
		 j < maxj_first; j++, aij += incaij, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = x_i[xi];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	    }
	    a_elem[0] = a_i[aij];
	    x_elem = x_i[xi];
	    prod[0] = x_elem * a_elem[0];
	    prod[1] = 0.0;
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	    aij += incaij2;
	    xi += incx;
	    for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

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
	      tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	      tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	    }

	    tmp1[0] = sum[0];
	    tmp1[1] = sum[1];
	    tmp1[0] = tmp2[0] + tmp1[0];
	    tmp1[1] = tmp2[1] + tmp1[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    if (i + 1 >= (n_i - k)) {
	      maxj_second--;
	    }
	    if (i >= k) {
	      astarti += (incaij + incaij2);
	      x_starti += incx;
	    } else {
	      maxj_first++;
	      astarti += incaij2;
	    }
	  }
	}
      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum[0] = sum[1] = 0.0;

	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_elem[1] = -a_elem[1];
	    x_elem = x_i[xi];
	    {
	      prod[0] = a_elem[0] * x_elem;
	      prod[1] = a_elem[1] * x_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  a_elem[0] = a_i[aij];
	  x_elem = x_i[xi];
	  prod[0] = x_elem * a_elem[0];
	  prod[1] = 0.0;
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	  aij += incaij2;
	  xi += incx;
	  for (j = 1; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];

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
	    tmp2[0] = y_elem[0] * beta_i[0] - y_elem[1] * beta_i[1];
	    tmp2[1] = y_elem[0] * beta_i[1] + y_elem[1] * beta_i[0];
	  }

	  {
	    tmp1[0] = sum[0] * alpha_i[0] - sum[1] * alpha_i[1];
	    tmp1[1] = sum[0] * alpha_i[1] + sum[1] * alpha_i[0];
	  }

	  tmp1[0] = tmp2[0] + tmp1[0];
	  tmp1[1] = tmp2[1] + tmp1[1];
	  y_i[yi] = tmp1[0];
	  y_i[yi + 1] = tmp1[1];
	  if (i + 1 >= (n_i - k)) {
	    maxj_second--;
	  }
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      }
    }
  }


}				/* end BLAS_chbmv_c_s */
