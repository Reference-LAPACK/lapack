#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dgemv_s_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta, double *y,
		    int incy)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * y, where A is a general matrix.
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
 * alpha        (input) double
 *              
 * A            (input) const float*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) double
 *
 * y            (input/output) double*
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
 */
{
  static const char routine_name[] = "BLAS_dgemv_s_d";

  int i, j;
  int iy, jx, kx, ky;
  int lenx, leny;
  int ai, aij;
  int incai, incaij;

  const float *a_i = a;
  const double *x_i = x;
  double *y_i = y;
  double alpha_i = alpha;
  double beta_i = beta;
  float a_elem;
  double x_elem;
  double y_elem;
  double prod;
  double sum;
  double tmp1;
  double tmp2;


  /* all error calls */
  if (m < 0)
    BLAS_error(routine_name, -3, m, 0);
  else if (n <= 0)
    BLAS_error(routine_name, -4, n, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -9, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -12, incy, 0);

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
  if ((order == blas_colmajor && lda < m) ||
      (order == blas_rowmajor && lda < n))
    BLAS_error(routine_name, -7, lda, NULL);








  if (incx > 0)
    kx = 0;
  else
    kx = (1 - lenx) * incx;
  if (incy > 0)
    ky = 0;
  else
    ky = (1 - leny) * incy;

  /* No extra-precision needed for alpha = 0 */
  if (alpha_i == 0.0) {
    if (beta_i == 0.0) {
      iy = ky;
      for (i = 0; i < leny; i++) {
	y_i[iy] = 0.0;
	iy += incy;
      }
    } else if (!(beta_i == 0.0)) {
      iy = ky;
      for (i = 0; i < leny; i++) {
	y_elem = y_i[iy];
	tmp1 = y_elem * beta_i;
	y_i[iy] = tmp1;
	iy += incy;
      }
    }
  } else {

    /* if beta = 0, we can save m multiplies: y = alpha*A*x */
    if (beta_i == 0.0) {
      /* save m more multiplies if alpha = 1 */
      if (alpha_i == 1.0) {
	ai = 0;
	iy = ky;
	for (i = 0; i < leny; i++) {
	  sum = 0.0;
	  aij = ai;
	  jx = kx;
	  for (j = 0; j < lenx; j++) {
	    a_elem = a_i[aij];

	    x_elem = x_i[jx];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	    aij += incaij;
	    jx += incx;
	  }
	  y_i[iy] = sum;
	  ai += incai;
	  iy += incy;
	}
      } else {
	ai = 0;
	iy = ky;
	for (i = 0; i < leny; i++) {
	  sum = 0.0;
	  aij = ai;
	  jx = kx;
	  for (j = 0; j < lenx; j++) {
	    a_elem = a_i[aij];

	    x_elem = x_i[jx];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	    aij += incaij;
	    jx += incx;
	  }
	  tmp1 = sum * alpha_i;
	  y_i[iy] = tmp1;
	  ai += incai;
	  iy += incy;
	}
      }
    } else {
      /* the most general form, y = alpha*A*x + beta*y */
      ai = 0;
      iy = ky;
      for (i = 0; i < leny; i++) {
	sum = 0.0;;
	aij = ai;
	jx = kx;
	for (j = 0; j < lenx; j++) {
	  a_elem = a_i[aij];

	  x_elem = x_i[jx];
	  prod = a_elem * x_elem;
	  sum = sum + prod;
	  aij += incaij;
	  jx += incx;
	}
	tmp1 = sum * alpha_i;
	y_elem = y_i[iy];
	tmp2 = y_elem * beta_i;
	tmp1 = tmp1 + tmp2;
	y_i[iy] = tmp1;
	ai += incai;
	iy += incy;
      }
    }

  }



}
