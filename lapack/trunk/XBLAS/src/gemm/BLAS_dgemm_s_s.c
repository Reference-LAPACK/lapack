#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dgemm_s_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    double alpha, const float *a, int lda, const float *b,
		    int ldb, double beta, double *c, int ldc)

/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) double
 *
 * a       (input) const float*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const float*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) double
 *
 * c       (input/output) double*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{
  static const char routine_name[] = "BLAS_dgemm_s_s";


  /* Integer Index Variables */
  int i, j, h;

  int ai, bj, ci;
  int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

  int incai, incaih;		/* Index increments for matrix a */
  int incbj, incbhj;		/* Index increments for matrix b */
  int incci, inccij;		/* Index increments for matrix c */

  /* Input Matrices */
  const float *a_i = a;
  const float *b_i = b;

  /* Output Matrix */
  double *c_i = c;

  /* Input Scalars */
  double alpha_i = alpha;
  double beta_i = beta;

  /* Temporary Floating-Point Variables */
  float a_elem;
  float b_elem;
  double c_elem;
  double prod;
  double sum;
  double tmp1;
  double tmp2;



  /* Test for error conditions */
  if (m < 0)
    BLAS_error(routine_name, -4, m, NULL);
  if (n < 0)
    BLAS_error(routine_name, -5, n, NULL);
  if (k < 0)
    BLAS_error(routine_name, -6, k, NULL);

  if (order == blas_colmajor) {

    if (ldc < m)
      BLAS_error(routine_name, -14, ldc, NULL);

    if (transa == blas_no_trans) {
      if (lda < m)
	BLAS_error(routine_name, -9, lda, NULL);
    } else {
      if (lda < k)
	BLAS_error(routine_name, -9, lda, NULL);
    }

    if (transb == blas_no_trans) {
      if (ldb < k)
	BLAS_error(routine_name, -11, ldb, NULL);
    } else {
      if (ldb < n)
	BLAS_error(routine_name, -11, ldb, NULL);
    }

  } else {
    /* row major */
    if (ldc < n)
      BLAS_error(routine_name, -14, ldc, NULL);

    if (transa == blas_no_trans) {
      if (lda < k)
	BLAS_error(routine_name, -9, lda, NULL);
    } else {
      if (lda < m)
	BLAS_error(routine_name, -9, lda, NULL);
    }

    if (transb == blas_no_trans) {
      if (ldb < n)
	BLAS_error(routine_name, -11, ldb, NULL);
    } else {
      if (ldb < k)
	BLAS_error(routine_name, -11, ldb, NULL);
    }
  }

  /* Test for no-op */
  if (n == 0 || m == 0 || k == 0)
    return;
  if (alpha_i == 0.0 && beta_i == 1.0) {
    return;
  }

  /* Set Index Parameters */
  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;

    if (transa == blas_no_trans) {
      incai = 1;
      incaih = lda;
    } else {
      incai = lda;
      incaih = 1;
    }

    if (transb == blas_no_trans) {
      incbj = ldb;
      incbhj = 1;
    } else {
      incbj = 1;
      incbhj = ldb;
    }

  } else {
    /* row major */
    incci = ldc;
    inccij = 1;

    if (transa == blas_no_trans) {
      incai = lda;
      incaih = 1;
    } else {
      incai = 1;
      incaih = lda;
    }

    if (transb == blas_no_trans) {
      incbj = 1;
      incbhj = ldb;
    } else {
      incbj = ldb;
      incbhj = 1;
    }

  }



  /* Ajustment to increments */







  /* alpha = 0.  In this case, just return beta * C */
  if (alpha_i == 0.0) {

    ci = 0;
    for (i = 0; i < m; i++, ci += incci) {
      cij = ci;
      for (j = 0; j < n; j++, cij += inccij) {
	c_elem = c_i[cij];
	tmp1 = c_elem * beta_i;
	c_i[cij] = tmp1;
      }
    }

  } else if (alpha_i == 1.0) {

    /* Case alpha == 1. */

    if (beta_i == 0.0) {
      /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

      ci = 0;
      ai = 0;
      for (i = 0; i < m; i++, ci += incci, ai += incai) {

	cij = ci;
	bj = 0;

	for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

	  aih = ai;
	  bhj = bj;

	  sum = 0.0;

	  for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	    a_elem = a_i[aih];
	    b_elem = b_i[bhj];
	    if (transa == blas_conj_trans) {

	    }
	    if (transb == blas_conj_trans) {

	    }
	    prod = (double) a_elem *b_elem;
	    sum = sum + prod;
	  }
	  c_i[cij] = sum;
	}
      }

    } else {
      /* Case alpha == 1, but beta != 0.
         We compute   C <--- A * B + beta * C   */

      ci = 0;
      ai = 0;
      for (i = 0; i < m; i++, ci += incci, ai += incai) {

	cij = ci;
	bj = 0;

	for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

	  aih = ai;
	  bhj = bj;

	  sum = 0.0;

	  for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	    a_elem = a_i[aih];
	    b_elem = b_i[bhj];
	    if (transa == blas_conj_trans) {

	    }
	    if (transb == blas_conj_trans) {

	    }
	    prod = (double) a_elem *b_elem;
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

    /* The most general form,   C <-- alpha * A * B + beta * C  */
    ci = 0;
    ai = 0;
    for (i = 0; i < m; i++, ci += incci, ai += incai) {

      cij = ci;
      bj = 0;

      for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

	aih = ai;
	bhj = bj;

	sum = 0.0;

	for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	  a_elem = a_i[aih];
	  b_elem = b_i[bhj];
	  if (transa == blas_conj_trans) {

	  }
	  if (transb == blas_conj_trans) {

	  }
	  prod = (double) a_elem *b_elem;
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
