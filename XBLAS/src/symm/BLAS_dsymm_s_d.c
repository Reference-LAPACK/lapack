#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dsymm_s_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const float *a, int lda,
		    const double *b, int ldb, double beta, double *c, int ldc)

/* 
 * Purpose
 * =======
 *
 * This routines computes one of the matrix product:
 *
 *     C  <-  alpha * A * B  +  beta * C
 *     C  <-  alpha * B * A  +  beta * C
 * 
 * where A is a symmetric matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input matrices A, B, and C.
 * 
 * side    (input) enum blas_side_type
 *         Determines which side of matrix B is matrix A is multiplied.
 *
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * m n     (input) int
 *         Size of matrices A, B, and C.
 *         Matrix A is m-by-m if it is multiplied on the left, 
 *                     n-by-n otherwise.
 *         Matrices B and C are m-by-n.
 *
 * alpha   (input) double
 * 
 * a       (input) const float*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * b       (input) const double*
 *         Matrix B.
 *   
 * ldb     (input) int
 *         Leading dimension of matrix B.
 *
 * beta    (input) double
 * 
 * c       (input/output) double*
 *         Matrix C.
 * 
 * ldc     (input) int
 *         Leading dimension of matrix C.
 *
 */
{

  /* Integer Index Variables */
  int i, j, k;

  int ai, bj, ci;
  int aik, bkj, cij;

  int incai, incbj, incci;
  int incaik1, incaik2, incbkj, inccij;

  int m_i, n_i;

  /* Input Matrices */
  const float *a_i = a;
  const double *b_i = b;

  /* Output Matrix */
  double *c_i = c;

  /* Input Scalars */
  double alpha_i = alpha;
  double beta_i = beta;

  /* Temporary Floating-Point Variables */
  float a_elem;
  double b_elem;
  double c_elem;
  double prod;
  double sum;
  double tmp1;
  double tmp2;



  /* Check for error conditions. */
  if (m <= 0 || n <= 0) {
    return;
  }

  if (order == blas_colmajor && (ldb < m || ldc < m)) {
    return;
  }
  if (order == blas_rowmajor && (ldb < n || ldc < n)) {
    return;
  }
  if (side == blas_left_side && lda < m) {
    return;
  }
  if (side == blas_right_side && lda < n) {
    return;
  }

  /* Test for no-op */
  if (alpha_i == 0.0 && beta_i == 1.0) {
    return;
  }

  /* Set Index Parameters */
  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incbj = ldb;
    incbkj = 1;
    incci = 1;
    inccij = ldc;
  } else {
    incbj = 1;
    incbkj = ldb;
    incci = ldc;
    inccij = 1;
  }

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaik1 = 1;
    incaik2 = lda;
  } else {
    incai = 1;
    incaik1 = lda;
    incaik2 = 1;
  }



  /* Adjustment to increments (if any) */








  /* alpha = 0.  In this case, just return beta * C */
  if (alpha_i == 0.0) {
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {
	c_elem = c_i[cij];
	tmp1 = c_elem * beta_i;
	c_i[cij] = tmp1;
      }
    }
  } else if (alpha_i == 1.0) {

    /* Case alpha == 1. */

    if (beta_i == 0.0) {
      /* Case alpha = 1, beta = 0.  We compute  C <--- A * B   or  B * A */
      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	for (j = 0, cij = ci, bj = 0; j < n_i;
	     j++, cij += inccij, bj += incbj) {

	  sum = 0.0;

	  for (k = 0, aik = ai, bkj = bj; k < i;
	       k++, aik += incaik1, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem = b_i[bkj];
	    prod = a_elem * b_elem;
	    sum = sum + prod;
	  }
	  for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem = b_i[bkj];
	    prod = a_elem * b_elem;
	    sum = sum + prod;
	  }
	  c_i[cij] = sum;
	}
      }
    } else {
      /* Case alpha = 1, but beta != 0. 
         We compute  C  <--- A * B + beta * C 
         or  C  <--- B * A + beta * C  */

      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	for (j = 0, cij = ci, bj = 0; j < n_i;
	     j++, cij += inccij, bj += incbj) {

	  sum = 0.0;

	  for (k = 0, aik = ai, bkj = bj; k < i;
	       k++, aik += incaik1, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem = b_i[bkj];
	    prod = a_elem * b_elem;
	    sum = sum + prod;
	  }
	  for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem = b_i[bkj];
	    prod = a_elem * b_elem;
	    sum = sum + prod;
	  }
	  c_elem = c_i[cij];
	  tmp2 = c_elem * beta_i;
	  tmp1 = sum;
	  tmp1 = tmp2 + tmp1;
	  c_i[cij] = tmp1;
	}
      }
    }

  } else {
    /* The most general form,   C <--- alpha * A * B + beta * C  
       or   C <--- alpha * B * A + beta * C  */

    for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
      for (j = 0, cij = ci, bj = 0; j < n_i; j++, cij += inccij, bj += incbj) {

	sum = 0.0;

	for (k = 0, aik = ai, bkj = bj; k < i;
	     k++, aik += incaik1, bkj += incbkj) {
	  a_elem = a_i[aik];
	  b_elem = b_i[bkj];
	  prod = a_elem * b_elem;
	  sum = sum + prod;
	}
	for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	  a_elem = a_i[aik];
	  b_elem = b_i[bkj];
	  prod = a_elem * b_elem;
	  sum = sum + prod;
	}
	tmp1 = sum * alpha_i;
	c_elem = c_i[cij];
	tmp2 = c_elem * beta_i;
	tmp1 = tmp1 + tmp2;
	c_i[cij] = tmp1;
      }
    }
  }



}
