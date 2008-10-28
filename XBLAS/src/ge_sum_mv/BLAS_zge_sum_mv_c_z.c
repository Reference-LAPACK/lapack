#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zge_sum_mv_c_z(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const void *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * B * y, 
 *     where A, B are general matricies.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * m            (input) int
 *              Row Dimension of A, B, length of output vector y
 *
 * n            (input) int
 *              Column Dimension of A, B and the length of vector x
 *
 * alpha        (input) const void*
 *              
 * A            (input) const void*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const void*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) const void*
 *
 * b            (input) const void*
 *
 * ldb          (input) int 
 *              Leading dimension of B
 *
 * y            (input/output) void*
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_zge_sum_mv_c_z";
  int i, j;
  int xi, yi;
  int x_starti, y_starti, incxi, incyi;
  int lda_min;
  int ai;
  int incai;
  int aij;
  int incaij;
  int bi;
  int incbi;
  int bij;
  int incbij;

  const float *a_i = (float *) a;
  const float *b_i = (float *) b;
  const double *x_i = (double *) x;
  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float a_elem[2];
  float b_elem[2];
  double x_elem[2];
  double prod[2];
  double sumA[2];
  double sumB[2];
  double tmp1[2];
  double tmp2[2];



  /* m is number of rows */
  /* n is number of columns */

  if (m == 0 || n == 0)
    return;


  /* all error calls */
  if (order == blas_rowmajor) {
    lda_min = n;
    incai = lda;		/* row stride */
    incbi = ldb;
    incbij = incaij = 1;	/* column stride */
  } else if (order == blas_colmajor) {
    lda_min = m;
    incai = incbi = 1;		/*row stride */
    incaij = lda;		/* column stride */
    incbij = ldb;
  } else {
    /* error, order not blas_colmajor not blas_rowmajor */
    BLAS_error(routine_name, -1, order, 0);
    return;
  }

  if (m < 0)
    BLAS_error(routine_name, -2, m, 0);
  else if (n < 0)
    BLAS_error(routine_name, -3, n, 0);
  if (lda < lda_min)
    BLAS_error(routine_name, -6, lda, 0);
  else if (ldb < lda_min)
    BLAS_error(routine_name, -11, ldb, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -8, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -13, incy, 0);

  incxi = incx;
  incyi = incy;
  incxi *= 2;
  incyi *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  if (incxi > 0)
    x_starti = 0;
  else
    x_starti = (1 - n) * incxi;

  if (incyi > 0)
    y_starti = 0;
  else
    y_starti = (1 - m) * incyi;



  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
      /* alpha, beta are 0.0 */
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	y_i[yi] = 0.0;
	y_i[yi + 1] = 0.0;
      }
    } else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
      /* alpha is 0.0, beta is 1.0 */


      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];

	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	y_i[yi] = sumB[0];
	y_i[yi + 1] = sumB[1];

	bi += incbi;
      }
    } else {
      /* alpha is 0.0, beta not 1.0 nor 0.0 */


      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];

	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	{
	  tmp1[0] =
	    (double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	  tmp1[1] =
	    (double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	}
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];

	bi += incbi;
      }
    }
  } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
    if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
      /* alpha is 1.0, beta is 0.0 */

      ai = 0;

      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;

	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;

	}
	/* now put the result into y_i */
	y_i[yi] = sumA[0];
	y_i[yi + 1] = sumA[1];
	ai += incai;

      }
    } else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
      /* alpha is 1.0, beta is 1.0 */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;
	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;
	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1[0] = sumA[0];
	tmp1[1] = sumA[1];
	tmp2[0] = sumB[0];
	tmp2[1] = sumB[1];
	tmp1[0] = tmp1[0] + tmp2[0];
	tmp1[1] = tmp1[1] + tmp2[1];
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
	ai += incai;
	bi += incbi;
      }
    } else {
      /* alpha is 1.0, beta is other */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;
	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;
	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1[0] = sumA[0];
	tmp1[1] = sumA[1];
	{
	  tmp2[0] =
	    (double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	  tmp2[1] =
	    (double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	}
	tmp1[0] = tmp1[0] + tmp2[0];
	tmp1[1] = tmp1[1] + tmp2[1];
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
	ai += incai;
	bi += incbi;
      }
    }
  } else {
    if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
      /* alpha is other, beta is 0.0 */

      ai = 0;

      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;

	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;

	}
	/* now put the result into y_i */
	{
	  tmp1[0] =
	    (double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	  tmp1[1] =
	    (double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	}
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
	ai += incai;

      }
    } else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
      /* alpha is other, beta is 1.0 */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;
	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;
	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	{
	  tmp1[0] =
	    (double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	  tmp1[1] =
	    (double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	}
	tmp2[0] = sumB[0];
	tmp2[1] = sumB[1];
	tmp1[0] = tmp1[0] + tmp2[0];
	tmp1[1] = tmp1[1] + tmp2[1];
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
	ai += incai;
	bi += incbi;
      }
    } else {
      /* most general form, alpha, beta are other */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA[0] = sumA[1] = 0.0;
	aij = ai;
	sumB[0] = sumB[1] = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  {
	    prod[0] =
	      (double) a_elem[0] * x_elem[0] - (double) a_elem[1] * x_elem[1];
	    prod[1] =
	      (double) a_elem[0] * x_elem[1] + (double) a_elem[1] * x_elem[0];
	  }
	  sumA[0] = sumA[0] + prod[0];
	  sumA[1] = sumA[1] + prod[1];
	  aij += incaij;
	  b_elem[0] = b_i[bij];
	  b_elem[1] = b_i[bij + 1];
	  {
	    prod[0] =
	      (double) b_elem[0] * x_elem[0] - (double) b_elem[1] * x_elem[1];
	    prod[1] =
	      (double) b_elem[0] * x_elem[1] + (double) b_elem[1] * x_elem[0];
	  }
	  sumB[0] = sumB[0] + prod[0];
	  sumB[1] = sumB[1] + prod[1];
	  bij += incbij;
	}
	/* now put the result into y_i */
	{
	  tmp1[0] =
	    (double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	  tmp1[1] =
	    (double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	}
	{
	  tmp2[0] =
	    (double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	  tmp2[1] =
	    (double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	}
	tmp1[0] = tmp1[0] + tmp2[0];
	tmp1[1] = tmp1[1] + tmp2[1];
	y_i[yi] = tmp1[0];
	y_i[yi + 1] = tmp1[1];
	ai += incai;
	bi += incbi;
      }
    }
  }


}				/* end BLAS_zge_sum_mv_c_z */
