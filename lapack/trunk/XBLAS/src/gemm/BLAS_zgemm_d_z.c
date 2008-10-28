#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zgemm_d_z(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const double *a, int lda,
		    const void *b, int ldb, const void *beta, void *c,
		    int ldc)

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
 * alpha   (input) const void*
 *
 * a       (input) const double*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const void*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) const void*
 *
 * c       (input/output) void*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{
  static const char routine_name[] = "BLAS_zgemm_d_z";


  /* Integer Index Variables */
  int i, j, h;

  int ai, bj, ci;
  int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

  int incai, incaih;		/* Index increments for matrix a */
  int incbj, incbhj;		/* Index increments for matrix b */
  int incci, inccij;		/* Index increments for matrix c */

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
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
      && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
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
  incci *= 2;
  inccij *= 2;


  incbj *= 2;
  incbhj *= 2;

  /* alpha = 0.  In this case, just return beta * C */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {

    ci = 0;
    for (i = 0; i < m; i++, ci += incci) {
      cij = ci;
      for (j = 0; j < n; j++, cij += inccij) {
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
      /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

      ci = 0;
      ai = 0;
      for (i = 0; i < m; i++, ci += incci, ai += incai) {

	cij = ci;
	bj = 0;

	for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

	  aih = ai;
	  bhj = bj;

	  sum[0] = sum[1] = 0.0;

	  for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	    a_elem = a_i[aih];
	    b_elem[0] = b_i[bhj];
	    b_elem[1] = b_i[bhj + 1];
	    if (transa == blas_conj_trans) {

	    }
	    if (transb == blas_conj_trans) {
	      b_elem[1] = -b_elem[1];
	    }
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

	  sum[0] = sum[1] = 0.0;

	  for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	    a_elem = a_i[aih];
	    b_elem[0] = b_i[bhj];
	    b_elem[1] = b_i[bhj + 1];
	    if (transa == blas_conj_trans) {

	    }
	    if (transb == blas_conj_trans) {
	      b_elem[1] = -b_elem[1];
	    }
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

    /* The most general form,   C <-- alpha * A * B + beta * C  */
    ci = 0;
    ai = 0;
    for (i = 0; i < m; i++, ci += incci, ai += incai) {

      cij = ci;
      bj = 0;

      for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

	aih = ai;
	bhj = bj;

	sum[0] = sum[1] = 0.0;

	for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	  a_elem = a_i[aih];
	  b_elem[0] = b_i[bhj];
	  b_elem[1] = b_i[bhj + 1];
	  if (transa == blas_conj_trans) {

	  }
	  if (transb == blas_conj_trans) {
	    b_elem[1] = -b_elem[1];
	  }
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
