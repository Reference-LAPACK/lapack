#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_csbmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const float *a,
		      int lda, const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a symmetric band matrix.
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
 * k       (input) int
 *         Number of subdiagonals ( = number of superdiagonals)
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
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 *
 *  Notes on storing a symmetric band matrix:
 * 
 *      Integers in the below arrays represent values of
 *              type float.
 *
 *    if we have a symettric matrix:
 *
 *      1  2  3  0  0
 *      2  4  5  6  0
 *      3  5  7  8  9
 *      0  6  8  10 11
 *      0  0  9  11 12
 *
 *     This matrix has n == 5, and k == 2. It can be stored in the
 *      following ways:
 *
 *      Notes for the examples:
 *      Each column below represents a contigous vector.
 *      Columns are strided by lda.
 *      An asterisk (*) represents a position in the 
 *       matrix that is not used.
 *      Note that the minimum lda (size of column) is 3 (k+1).
 *       lda may be arbitrarily large; an lda > 3 would mean
 *       there would be unused data at the bottom of the below
 *       columns.        
 *
 *    blas_colmajor and blas_upper:
 *      *  *  3  6  9  
 *      *  2  5  8  11 
 *      1  4  7  10 12
 *
 *
 *    blas_colmajor and blas_lower
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_upper 
 *      Columns here also represent contigous arrays.
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_lower
 *      Columns here also represent contigous arrays.
 *      *  *  3  6   9
 *      *  2  5  8   11
 *      1  4  7  10  12
 *
 */
{
  static const char routine_name[] = "BLAS_csbmv_s_s_x";
  switch (prec) {

  case blas_prec_single:{

      /* Integer Index Variables */
      int i, j;
      int xi, yi;
      int aij, astarti, x_starti, y_starti;
      int incaij, incaij2;
      int n_i;
      int maxj_first, maxj_second;

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
	incaij = 1;		/* increment in first loop */
	incaij2 = lda - 1;	/* increment in second loop */
	astarti = k;		/* does not start on zero element */
      } else {
	incaij = lda - 1;
	incaij2 = 1;
	astarti = 0;		/* start on first element of array */
      }
      /* Adjustment to increments (if any) */
      incy *= 2;
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
	/* number of elements done in first loop 
	   (this will increase by one over each column up to a limit) */
	maxj_first = 0;
	/* number of elements done in second loop the first time */
	maxj_second = MIN(k + 1, n_i);

	/* Case alpha == 1. */
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;
	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
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
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha != 1, but beta == 0. 
	       We compute  y  <--- A * x * a */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      y_elem[0] = y_i[yi];
	      y_elem[1] = y_i[yi + 1];
	      {
		tmp1[0] = alpha_i[0] * sum;
		tmp1[1] = alpha_i[1] * sum;
	      }
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
	      if (i >= k) {
		astarti += (incaij + incaij2);
		x_starti += incx;
	      } else {
		maxj_first++;
		astarti += incaij2;
	      }
	    }
	  } else {
	    /* The most general form,   y <--- alpha * A * x + beta * y */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = a_elem * x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
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
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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


      break;
    }
  case blas_prec_indigenous:
  case blas_prec_double:{

      /* Integer Index Variables */
      int i, j;
      int xi, yi;
      int aij, astarti, x_starti, y_starti;
      int incaij, incaij2;
      int n_i;
      int maxj_first, maxj_second;

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
      double prod;
      double sum;
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
	incaij = 1;		/* increment in first loop */
	incaij2 = lda - 1;	/* increment in second loop */
	astarti = k;		/* does not start on zero element */
      } else {
	incaij = lda - 1;
	incaij2 = 1;
	astarti = 0;		/* start on first element of array */
      }
      /* Adjustment to increments (if any) */
      incy *= 2;
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
	/*  determine the loop interation counts */
	/* number of elements done in first loop 
	   (this will increase by one over each column up to a limit) */
	maxj_first = 0;
	/* number of elements done in second loop the first time */
	maxj_second = MIN(k + 1, n_i);

	/* Case alpha == 1. */
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;
	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
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
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      tmp1[0] = tmp2[0] + tmp1[0];
	      tmp1[1] = tmp2[1] + tmp1[1];
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha != 1, but beta == 0. 
	       We compute  y  <--- A * x * a */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      y_elem[0] = y_i[yi];
	      y_elem[1] = y_i[yi + 1];
	      {
		tmp1[0] = alpha_i[0] * sum;
		tmp1[1] = alpha_i[1] * sum;
	      }
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
	      if (i >= k) {
		astarti += (incaij + incaij2);
		x_starti += incx;
	      } else {
		maxj_first++;
		astarti += incaij2;
	      }
	    }
	  } else {
	    /* The most general form,   y <--- alpha * A * x + beta * y */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		prod = (double) a_elem *x_elem;
		sum = sum + prod;
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
	      {
		tmp1[0] = alpha_i[0] * sum;
		tmp1[1] = alpha_i[1] * sum;
	      }
	      tmp1[0] = tmp2[0] + tmp1[0];
	      tmp1[1] = tmp2[1] + tmp1[1];
	      y_i[yi] = tmp1[0];
	      y_i[yi + 1] = tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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


      break;
    }

  case blas_prec_extra:{

      /* Integer Index Variables */
      int i, j;
      int xi, yi;
      int aij, astarti, x_starti, y_starti;
      int incaij, incaij2;
      int n_i;
      int maxj_first, maxj_second;

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
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
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
	incaij = 1;		/* increment in first loop */
	incaij2 = lda - 1;	/* increment in second loop */
	astarti = k;		/* does not start on zero element */
      } else {
	incaij = lda - 1;
	incaij2 = 1;
	astarti = 0;		/* start on first element of array */
      }
      /* Adjustment to increments (if any) */
      incy *= 2;
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

      FPU_FIX_START;

      /* alpha = 0.  In this case, just return beta * y */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
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
	    head_tmp1[0] = head_e1;
	    tail_tmp1[0] = tail_e1;
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
	    head_tmp1[1] = head_e1;
	    tail_tmp1[1] = tail_e1;
	  }
	  y_i[yi] = head_tmp1[0];
	  y_i[yi + 1] = head_tmp1[1];
	}
      } else {
	/*  determine the loop interation counts */
	/* number of elements done in first loop 
	   (this will increase by one over each column up to a limit) */
	maxj_first = 0;
	/* number of elements done in second loop the first time */
	maxj_second = MIN(k + 1, n_i);

	/* Case alpha == 1. */
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      head_sum = tail_sum = 0.0;
	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      head_tmp1[0] = head_sum;
	      tail_tmp1[0] = tail_sum;
	      head_tmp1[1] = tail_tmp1[1] = 0.0;
	      y_i[yi] = head_tmp1[0];
	      y_i[yi + 1] = head_tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	      head_sum = tail_sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      head_tmp1[0] = head_sum;
	      tail_tmp1[0] = tail_sum;
	      head_tmp1[1] = tail_tmp1[1] = 0.0;
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_tmp2[0];
		tail_a = tail_tmp2[0];
		head_b = head_tmp1[0];
		tail_b = tail_tmp1[0];
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
		head_a = head_tmp2[1];
		tail_a = tail_tmp2[1];
		head_b = head_tmp1[1];
		tail_b = tail_tmp1[1];
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
	      y_i[yi] = head_tmp1[0];
	      y_i[yi + 1] = head_tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    /* Case alpha != 1, but beta == 0. 
	       We compute  y  <--- A * x * a */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      head_sum = tail_sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      y_elem[0] = y_i[yi];
	      y_elem[1] = y_i[yi + 1];
	      {
		double head_e1, tail_e1;
		double dt;
		dt = (double) alpha_i[0];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_sum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_e1 = t1 + t2;
		  tail_e1 = t2 - (head_e1 - t1);
		}
		head_tmp1[0] = head_e1;
		tail_tmp1[0] = tail_e1;
		dt = (double) alpha_i[1];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_sum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_e1 = t1 + t2;
		  tail_e1 = t2 - (head_e1 - t1);
		}
		head_tmp1[1] = head_e1;
		tail_tmp1[1] = tail_e1;
	      }
	      y_i[yi] = head_tmp1[0];
	      y_i[yi + 1] = head_tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
	      if (i >= k) {
		astarti += (incaij + incaij2);
		x_starti += incx;
	      } else {
		maxj_first++;
		astarti += incaij2;
	      }
	    }
	  } else {
	    /* The most general form,   y <--- alpha * A * x + beta * y */
	    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	      head_sum = tail_sum = 0.0;

	      for (j = 0, aij = astarti, xi = x_starti;
		   j < maxj_first; j++, aij += incaij, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
	      }
	      for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
		a_elem = a_i[aij];
		x_elem = x_i[xi];
		head_prod = (double) a_elem *x_elem;
		tail_prod = 0.0;
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
		double head_e1, tail_e1;
		double dt;
		dt = (double) alpha_i[0];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_sum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_e1 = t1 + t2;
		  tail_e1 = t2 - (head_e1 - t1);
		}
		head_tmp1[0] = head_e1;
		tail_tmp1[0] = tail_e1;
		dt = (double) alpha_i[1];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_sum * split;
		  a11 = con - head_sum;
		  a11 = con - a11;
		  a21 = head_sum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_sum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_sum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_e1 = t1 + t2;
		  tail_e1 = t2 - (head_e1 - t1);
		}
		head_tmp1[1] = head_e1;
		tail_tmp1[1] = tail_e1;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_tmp2[0];
		tail_a = tail_tmp2[0];
		head_b = head_tmp1[0];
		tail_b = tail_tmp1[0];
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
		head_a = head_tmp2[1];
		tail_a = tail_tmp2[1];
		head_b = head_tmp1[1];
		tail_b = tail_tmp1[1];
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
	      y_i[yi] = head_tmp1[0];
	      y_i[yi + 1] = head_tmp1[1];
	      if (i + 1 >= (n_i - k))
		maxj_second--;
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
      FPU_FIX_STOP;

      break;
    }

  default:
    BLAS_error(routine_name, -13, prec, 0);
    break;
  }
}				/* end BLAS_csbmv_s_s */
