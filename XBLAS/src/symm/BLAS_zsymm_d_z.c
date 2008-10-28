#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zsymm_d_z(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const double *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc)

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
 * alpha   (input) const void*
 * 
 * a       (input) const double*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * b       (input) const void*
 *         Matrix B.
 *   
 * ldb     (input) int
 *         Leading dimension of matrix B.
 *
 * beta    (input) const void*
 * 
 * c       (input/output) void*
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
  const double *a_i = a;
  const double *b_i = (double *) b;

  /* Output Matrix */
  double *c_i = (double *) c;

  /* Input Scalars */
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary Floating-Point Variables */
  double a_elem;
  double b_elem[2];
  double c_elem[2];
  double prod[2];
  double sum[2];
  double tmp1[2];
  double tmp2[2];



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
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
      && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
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
  incci *= 2;
  inccij *= 2;



  incbj *= 2;
  incbkj *= 2;

  /* alpha = 0.  In this case, just return beta * C */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[cij];
	c_elem[1] = c_i[cij + 1];
	{
	  tmp1[0] =
	    (double) c_elem[0] * beta_i[0] - (double) c_elem[1] * beta_i[1];
	  tmp1[1] =
	    (double) c_elem[0] * beta_i[1] + (double) c_elem[1] * beta_i[0];
	}
	c_i[cij] = tmp1[0];
	c_i[cij + 1] = tmp1[1];
      }
    }
  } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

    /* Case alpha == 1. */

    if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
      /* Case alpha = 1, beta = 0.  We compute  C <--- A * B   or  B * A */
      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	for (j = 0, cij = ci, bj = 0; j < n_i;
	     j++, cij += inccij, bj += incbj) {

	  sum[0] = sum[1] = 0.0;

	  for (k = 0, aik = ai, bkj = bj; k < i;
	       k++, aik += incaik1, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem[0] = b_i[bkj];
	    b_elem[1] = b_i[bkj + 1];
	    {
	      prod[0] = b_elem[0] * a_elem;
	      prod[1] = b_elem[1] * a_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem[0] = b_i[bkj];
	    b_elem[1] = b_i[bkj + 1];
	    {
	      prod[0] = b_elem[0] * a_elem;
	      prod[1] = b_elem[1] * a_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  c_i[cij] = sum[0];
	  c_i[cij + 1] = sum[1];
	}
      }
    } else {
      /* Case alpha = 1, but beta != 0. 
         We compute  C  <--- A * B + beta * C 
         or  C  <--- B * A + beta * C  */

      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	for (j = 0, cij = ci, bj = 0; j < n_i;
	     j++, cij += inccij, bj += incbj) {

	  sum[0] = sum[1] = 0.0;

	  for (k = 0, aik = ai, bkj = bj; k < i;
	       k++, aik += incaik1, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem[0] = b_i[bkj];
	    b_elem[1] = b_i[bkj + 1];
	    {
	      prod[0] = b_elem[0] * a_elem;
	      prod[1] = b_elem[1] * a_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	    a_elem = a_i[aik];
	    b_elem[0] = b_i[bkj];
	    b_elem[1] = b_i[bkj + 1];
	    {
	      prod[0] = b_elem[0] * a_elem;
	      prod[1] = b_elem[1] * a_elem;
	    }
	    sum[0] = sum[0] + prod[0];
	    sum[1] = sum[1] + prod[1];
	  }
	  c_elem[0] = c_i[cij];
	  c_elem[1] = c_i[cij + 1];
	  {
	    tmp2[0] =
	      (double) c_elem[0] * beta_i[0] - (double) c_elem[1] * beta_i[1];
	    tmp2[1] =
	      (double) c_elem[0] * beta_i[1] + (double) c_elem[1] * beta_i[0];
	  }
	  tmp1[0] = sum[0];
	  tmp1[1] = sum[1];
	  tmp1[0] = tmp2[0] + tmp1[0];
	  tmp1[1] = tmp2[1] + tmp1[1];
	  c_i[cij] = tmp1[0];
	  c_i[cij + 1] = tmp1[1];
	}
      }
    }

  } else {
    /* The most general form,   C <--- alpha * A * B + beta * C  
       or   C <--- alpha * B * A + beta * C  */

    for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
      for (j = 0, cij = ci, bj = 0; j < n_i; j++, cij += inccij, bj += incbj) {

	sum[0] = sum[1] = 0.0;

	for (k = 0, aik = ai, bkj = bj; k < i;
	     k++, aik += incaik1, bkj += incbkj) {
	  a_elem = a_i[aik];
	  b_elem[0] = b_i[bkj];
	  b_elem[1] = b_i[bkj + 1];
	  {
	    prod[0] = b_elem[0] * a_elem;
	    prod[1] = b_elem[1] * a_elem;
	  }
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	}
	for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	  a_elem = a_i[aik];
	  b_elem[0] = b_i[bkj];
	  b_elem[1] = b_i[bkj + 1];
	  {
	    prod[0] = b_elem[0] * a_elem;
	    prod[1] = b_elem[1] * a_elem;
	  }
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];
	}
	{
	  tmp1[0] =
	    (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	  tmp1[1] =
	    (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	}
	c_elem[0] = c_i[cij];
	c_elem[1] = c_i[cij + 1];
	{
	  tmp2[0] =
	    (double) c_elem[0] * beta_i[0] - (double) c_elem[1] * beta_i[1];
	  tmp2[1] =
	    (double) c_elem[0] * beta_i[1] + (double) c_elem[1] * beta_i[0];
	}
	tmp1[0] = tmp1[0] + tmp2[0];
	tmp1[1] = tmp1[1] + tmp2[1];
	c_i[cij] = tmp1[0];
	c_i[cij + 1] = tmp1[1];
      }
    }
  }



}
