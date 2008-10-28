#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_cgemm_s_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const float *a, int lda,
		      const float *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec)

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
 * beta    (input) const void*
 *
 * c       (input/output) void*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  static const char routine_name[] = "BLAS_cgemm_s_s_x";
  switch (prec) {

  case blas_prec_single:{


      /* Integer Index Variables */
      int i, j, h;

      int ai, bj, ci;
      int aih, bhj, cij;	/* Index into matrices a, b, c during multiply */

      int incai, incaih;	/* Index increments for matrix a */
      int incbj, incbhj;	/* Index increments for matrix b */
      int incci, inccij;	/* Index increments for matrix c */

      /* Input Matrices */
      const float *a_i = a;
      const float *b_i = b;

      /* Output Matrix */
      float *c_i = (float *) c;

      /* Input Scalars */
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;

      /* Temporary Floating-Point Variables */
      float a_elem;
      float b_elem;
      float c_elem[2];
      float prod;
      float sum;
      float tmp1[2];
      float tmp2[2];



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





      /* alpha = 0.  In this case, just return beta * C */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci) {
	  cij = ci;
	  for (j = 0; j < n; j++, cij += inccij) {
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      tmp1[0] = c_elem[0] * beta_i[0] - c_elem[1] * beta_i[1];
	      tmp1[1] = c_elem[0] * beta_i[1] + c_elem[1] * beta_i[0];
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

	      sum = 0.0;

	      for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
		a_elem = a_i[aih];
		b_elem = b_i[bhj];
		if (transa == blas_conj_trans) {

		}
		if (transb == blas_conj_trans) {

		}
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      c_i[cij] = tmp1[0];
	      c_i[cij + 1] = tmp1[1];
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
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }

	      c_elem[0] = c_i[cij];
	      c_elem[1] = c_i[cij + 1];
	      {
		tmp2[0] = c_elem[0] * beta_i[0] - c_elem[1] * beta_i[1];
		tmp2[1] = c_elem[0] * beta_i[1] + c_elem[1] * beta_i[0];
	      }

	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
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

	    sum = 0.0;

	    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	      a_elem = a_i[aih];
	      b_elem = b_i[bhj];
	      if (transa == blas_conj_trans) {

	      }
	      if (transb == blas_conj_trans) {

	      }
	      prod = a_elem * b_elem;
	      sum = sum + prod;
	    }

	    {
	      tmp1[0] = alpha_i[0] * sum;
	      tmp1[1] = alpha_i[1] * sum;
	    }
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      tmp2[0] = c_elem[0] * beta_i[0] - c_elem[1] * beta_i[1];
	      tmp2[1] = c_elem[0] * beta_i[1] + c_elem[1] * beta_i[0];
	    }

	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    c_i[cij] = tmp1[0];
	    c_i[cij + 1] = tmp1[1];
	  }
	}

      }



      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:{


      /* Integer Index Variables */
      int i, j, h;

      int ai, bj, ci;
      int aih, bhj, cij;	/* Index into matrices a, b, c during multiply */

      int incai, incaih;	/* Index increments for matrix a */
      int incbj, incbhj;	/* Index increments for matrix b */
      int incci, inccij;	/* Index increments for matrix c */

      /* Input Matrices */
      const float *a_i = a;
      const float *b_i = b;

      /* Output Matrix */
      float *c_i = (float *) c;

      /* Input Scalars */
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;

      /* Temporary Floating-Point Variables */
      float a_elem;
      float b_elem;
      float c_elem[2];
      double prod;
      double sum;
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
		(double) c_elem[0] * beta_i[0] -
		(double) c_elem[1] * beta_i[1];
	      tmp1[1] =
		(double) c_elem[0] * beta_i[1] +
		(double) c_elem[1] * beta_i[0];
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
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      c_i[cij] = tmp1[0];
	      c_i[cij + 1] = tmp1[1];
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

	      c_elem[0] = c_i[cij];
	      c_elem[1] = c_i[cij + 1];
	      {
		tmp2[0] =
		  (double) c_elem[0] * beta_i[0] -
		  (double) c_elem[1] * beta_i[1];
		tmp2[1] =
		  (double) c_elem[0] * beta_i[1] +
		  (double) c_elem[1] * beta_i[0];
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
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

	    {
	      tmp1[0] = alpha_i[0] * sum;
	      tmp1[1] = alpha_i[1] * sum;
	    }
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      tmp2[0] =
		(double) c_elem[0] * beta_i[0] -
		(double) c_elem[1] * beta_i[1];
	      tmp2[1] =
		(double) c_elem[0] * beta_i[1] +
		(double) c_elem[1] * beta_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    c_i[cij] = tmp1[0];
	    c_i[cij + 1] = tmp1[1];
	  }
	}

      }



      break;
    }

  case blas_prec_extra:{


      /* Integer Index Variables */
      int i, j, h;

      int ai, bj, ci;
      int aih, bhj, cij;	/* Index into matrices a, b, c during multiply */

      int incai, incaih;	/* Index increments for matrix a */
      int incbj, incbhj;	/* Index increments for matrix b */
      int incci, inccij;	/* Index increments for matrix c */

      /* Input Matrices */
      const float *a_i = a;
      const float *b_i = b;

      /* Output Matrix */
      float *c_i = (float *) c;

      /* Input Scalars */
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;

      /* Temporary Floating-Point Variables */
      float a_elem;
      float b_elem;
      float c_elem[2];
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];

      FPU_FIX_DECL;

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

      FPU_FIX_START;

      /* Ajustment to increments */
      incci *= 2;
      inccij *= 2;





      /* alpha = 0.  In this case, just return beta * C */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci) {
	  cij = ci;
	  for (j = 0; j < n; j++, cij += inccij) {
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) c_elem[0] * beta_i[0];
	      d2 = (double) -c_elem[1] * beta_i[1];
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
	      d1 = (double) c_elem[0] * beta_i[1];
	      d2 = (double) c_elem[1] * beta_i[0];
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
	    c_i[cij] = head_tmp1[0];
	    c_i[cij + 1] = head_tmp1[1];
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

	      head_sum = tail_sum = 0.0;

	      for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
		a_elem = a_i[aih];
		b_elem = b_i[bhj];
		if (transa == blas_conj_trans) {

		}
		if (transb == blas_conj_trans) {

		}
		head_prod = (double) a_elem *b_elem;
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
	      c_i[cij] = head_tmp1[0];
	      c_i[cij + 1] = head_tmp1[1];
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

	      head_sum = tail_sum = 0.0;

	      for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
		a_elem = a_i[aih];
		b_elem = b_i[bhj];
		if (transa == blas_conj_trans) {

		}
		if (transb == blas_conj_trans) {

		}
		head_prod = (double) a_elem *b_elem;
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

	      c_elem[0] = c_i[cij];
	      c_elem[1] = c_i[cij + 1];
	      {
		double head_e1, tail_e1;
		double d1;
		double d2;
		/* Real part */
		d1 = (double) c_elem[0] * beta_i[0];
		d2 = (double) -c_elem[1] * beta_i[1];
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
		d1 = (double) c_elem[0] * beta_i[1];
		d2 = (double) c_elem[1] * beta_i[0];
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
	      c_i[cij] = head_tmp1[0];
	      c_i[cij + 1] = head_tmp1[1];
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

	    head_sum = tail_sum = 0.0;

	    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
	      a_elem = a_i[aih];
	      b_elem = b_i[bhj];
	      if (transa == blas_conj_trans) {

	      }
	      if (transb == blas_conj_trans) {

	      }
	      head_prod = (double) a_elem *b_elem;
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
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      double head_e1, tail_e1;
	      double d1;
	      double d2;
	      /* Real part */
	      d1 = (double) c_elem[0] * beta_i[0];
	      d2 = (double) -c_elem[1] * beta_i[1];
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
	      d1 = (double) c_elem[0] * beta_i[1];
	      d2 = (double) c_elem[1] * beta_i[0];
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
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
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
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
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
	    c_i[cij] = head_tmp1[0];
	    c_i[cij + 1] = head_tmp1[1];
	  }
	}

      }

      FPU_FIX_STOP;

      break;
    }
  }
}
