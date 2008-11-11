#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_zsymv2_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const double *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     const void *beta, void *y, int incy)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a symmetric matrix.
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
 * a       (input) double*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) double*
 *         Vector x_head
 *
 * x_tail  (input) double*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) const void*
 * 
 * y       (input) double*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
 *
 */
{
  /* Routine name */
  const char routine_name[] = "BLAS_zsymv2_d_d";

  int i, j;
  int xi, yi, xi0, yi0;
  int aij, ai;
  int incai;
  int incaij, incaij2;

  const double *a_i = a;
  const double *x_head_i = x_head;
  const double *x_tail_i = x_tail;
  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double a_elem;
  double x_elem;
  double y_elem[2];
  double prod1;
  double prod2;
  double sum1;
  double sum2;
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


  incy *= 2;



  xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
  yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



  /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
  for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
    sum1 = 0.0;
    sum2 = 0.0;

    for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
      a_elem = a_i[aij];
      x_elem = x_head_i[xi];
      prod1 = a_elem * x_elem;
      sum1 = sum1 + prod1;
      x_elem = x_tail_i[xi];
      prod2 = a_elem * x_elem;
      sum2 = sum2 + prod2;
    }
    for (; j < n; j++, aij += incaij2, xi += incx) {
      a_elem = a_i[aij];
      x_elem = x_head_i[xi];
      prod1 = a_elem * x_elem;
      sum1 = sum1 + prod1;
      x_elem = x_tail_i[xi];
      prod2 = a_elem * x_elem;
      sum2 = sum2 + prod2;
    }
    sum1 = sum1 + sum2;
    {
      tmp1[0] = alpha_i[0] * sum1;
      tmp1[1] = alpha_i[1] * sum1;
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



}				/* end BLAS_zsymv2_d_d */
