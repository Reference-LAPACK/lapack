#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zgemv2_z_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const double *head_x, const double *tail_x, int incx,
		     const void *beta, void *y, int incy)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * op(A) * head_x + alpha * op(A) * tail_x + beta * y,
 * where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of A: no trans, trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of A
 *
 * n            (input) int
 *              Dimension of A and the length of vector x and z
 *
 * alpha        (input) const void*
 *              
 * A            (input) const void*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * head_x
 * tail_x       (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) const void*
 *
 * y            (input) const void*
 * 
 * incy         (input) int
 *              The stride for vector y.
 * 
 */
{
  static const char routine_name[] = "BLAS_zgemv2_z_d";

  int i, j;
  int iy, jx, kx, ky;
  int lenx, leny;
  int ai, aij;
  int incai, incaij;

  const double *a_i = (double *) a;
  const double *head_x_i = head_x;
  const double *tail_x_i = tail_x;
  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double a_elem[2];
  double x_elem;
  double y_elem[2];
  double prod[2];
  double sum[2];
  double sum2[2];
  double tmp1[2];
  double tmp2[2];


  /* all error calls */
  if (m < 0)
    BLAS_error(routine_name, -3, m, 0);
  else if (n <= 0)
    BLAS_error(routine_name, -4, n, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -10, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -13, incy, 0);

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

  if (lda < leny)
    BLAS_error(routine_name, -7, lda, NULL);




  incy *= 2;
  incai *= 2;
  incaij *= 2;

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
	    (double) y_elem[0] * beta_i[0] - (double) y_elem[1] * beta_i[1];
	  tmp1[1] =
	    (double) y_elem[0] * beta_i[1] + (double) y_elem[1] * beta_i[0];
	}
	y_i[iy] = tmp1[0];
	y_i[iy + 1] = tmp1[1];
	iy += incy;
      }
    }
  } else {			/* alpha != 0 */
    if (trans == blas_conj_trans) {

      /* if beta = 0, we can save m multiplies:
         y = alpha*A*head_x + alpha*A*tail_x  */
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  /* save m more multiplies if alpha = 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;
	    sum2[0] = sum2[1] = 0.0;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    sum[0] = sum[0] + sum2[0];
	    sum[1] = sum[1] + sum2[1];
	    y_i[iy] = sum[0];
	    y_i[iy + 1] = sum[1];
	    ai += incai;
	    iy += incy;
	  }			/* end for */
	} else {		/* alpha != 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;
	    sum2[0] = sum2[1] = 0.0;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      tmp1[0] =
		(double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	      tmp1[1] =
		(double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	    }
	    {
	      tmp2[0] =
		(double) sum2[0] * alpha_i[0] - (double) sum2[1] * alpha_i[1];
	      tmp2[1] =
		(double) sum2[0] * alpha_i[1] + (double) sum2[1] * alpha_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[iy] = tmp1[0];
	    y_i[iy + 1] = tmp1[1];
	    ai += incai;
	    iy += incy;
	  }
	}
      } else {			/* beta != 0 */
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  /* save m multiplies if alpha = 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;;
	    sum2[0] = sum2[1] = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    sum[0] = sum[0] + sum2[0];
	    sum[1] = sum[1] + sum2[1];
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
	    tmp2[0] = sum[0] + tmp1[0];
	    tmp2[1] = sum[1] + tmp1[1];
	    y_i[iy] = tmp2[0];
	    y_i[iy + 1] = tmp2[1];
	    ai += incai;
	    iy += incy;
	  }
	} else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;;
	    sum2[0] = sum2[1] = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      a_elem[1] = -a_elem[1];
	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      tmp1[0] =
		(double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	      tmp1[1] =
		(double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	    }
	    {
	      tmp2[0] =
		(double) sum2[0] * alpha_i[0] - (double) sum2[1] * alpha_i[1];
	      tmp2[1] =
		(double) sum2[0] * alpha_i[1] + (double) sum2[1] * alpha_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
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

    } else {

      /* if beta = 0, we can save m multiplies:
         y = alpha*A*head_x + alpha*A*tail_x  */
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  /* save m more multiplies if alpha = 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;
	    sum2[0] = sum2[1] = 0.0;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    sum[0] = sum[0] + sum2[0];
	    sum[1] = sum[1] + sum2[1];
	    y_i[iy] = sum[0];
	    y_i[iy + 1] = sum[1];
	    ai += incai;
	    iy += incy;
	  }			/* end for */
	} else {		/* alpha != 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;
	    sum2[0] = sum2[1] = 0.0;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      tmp1[0] =
		(double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	      tmp1[1] =
		(double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	    }
	    {
	      tmp2[0] =
		(double) sum2[0] * alpha_i[0] - (double) sum2[1] * alpha_i[1];
	      tmp2[1] =
		(double) sum2[0] * alpha_i[1] + (double) sum2[1] * alpha_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[iy] = tmp1[0];
	    y_i[iy + 1] = tmp1[1];
	    ai += incai;
	    iy += incy;
	  }
	}
      } else {			/* beta != 0 */
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  /* save m multiplies if alpha = 1 */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;;
	    sum2[0] = sum2[1] = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    sum[0] = sum[0] + sum2[0];
	    sum[1] = sum[1] + sum2[1];
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
	    tmp2[0] = sum[0] + tmp1[0];
	    tmp2[1] = sum[1] + tmp1[1];
	    y_i[iy] = tmp2[0];
	    y_i[iy + 1] = tmp2[1];
	    ai += incai;
	    iy += incy;
	  }
	} else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum[0] = sum[1] = 0.0;;
	    sum2[0] = sum2[1] = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];

	      x_elem = head_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];
	      x_elem = tail_x_i[jx];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sum2[0] = sum2[0] + prod[0];
	      sum2[1] = sum2[1] + prod[1];
	      aij += incaij;
	      jx += incx;
	    }
	    {
	      tmp1[0] =
		(double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	      tmp1[1] =
		(double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	    }
	    {
	      tmp2[0] =
		(double) sum2[0] * alpha_i[0] - (double) sum2[1] * alpha_i[1];
	      tmp2[1] =
		(double) sum2[0] * alpha_i[1] + (double) sum2[1] * alpha_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
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

    }
  }



}
