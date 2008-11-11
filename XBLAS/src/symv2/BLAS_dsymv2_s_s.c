#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_dsymv2_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const float *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     double beta, double *y, int incy)

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
 * alpha   (input) double
 * 
 * a       (input) float*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) float*
 *         Vector x_head
 *
 * x_tail  (input) float*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) double
 * 
 * y       (input) float*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
 *
 */
{
  /* Routine name */
  const char routine_name[] = "BLAS_dsymv2_s_s";

  int i, j;
  int xi, yi, xi0, yi0;
  int aij, ai;
  int incai;
  int incaij, incaij2;

  const float *a_i = a;
  const float *x_head_i = x_head;
  const float *x_tail_i = x_tail;
  double *y_i = y;
  double alpha_i = alpha;
  double beta_i = beta;
  float a_elem;
  float x_elem;
  double y_elem;
  double prod1;
  double prod2;
  double sum1;
  double sum2;
  double tmp1;
  double tmp2;
  double tmp3;



  /* Test for no-op */
  if (n <= 0) {
    return;
  }
  if (alpha_i == 0.0 && beta_i == 1.0) {
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






  xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
  yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



  /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
  for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
    sum1 = 0.0;
    sum2 = 0.0;

    for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
      a_elem = a_i[aij];
      x_elem = x_head_i[xi];
      prod1 = (double) a_elem *x_elem;
      sum1 = sum1 + prod1;
      x_elem = x_tail_i[xi];
      prod2 = (double) a_elem *x_elem;
      sum2 = sum2 + prod2;
    }
    for (; j < n; j++, aij += incaij2, xi += incx) {
      a_elem = a_i[aij];
      x_elem = x_head_i[xi];
      prod1 = (double) a_elem *x_elem;
      sum1 = sum1 + prod1;
      x_elem = x_tail_i[xi];
      prod2 = (double) a_elem *x_elem;
      sum2 = sum2 + prod2;
    }
    sum1 = sum1 + sum2;
    tmp1 = sum1 * alpha_i;
    y_elem = y_i[yi];
    tmp2 = y_elem * beta_i;
    tmp3 = tmp1 + tmp2;
    y_i[yi] = tmp3;
  }



}				/* end BLAS_dsymv2_s_s */
