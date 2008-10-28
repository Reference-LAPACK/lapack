#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dgemm_d_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      double alpha, const double *a, int lda, const float *b,
		      int ldb, double beta, double *c, int ldc,
		      enum blas_prec_type prec)

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
 * a       (input) const double*
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
  static const char routine_name[] = "BLAS_dgemm_d_s_x";
  switch (prec) {

  case blas_prec_single:
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
      const double *a_i = a;
      const float *b_i = b;

      /* Output Matrix */
      double *c_i = c;

      /* Input Scalars */
      double alpha_i = alpha;
      double beta_i = beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
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
		prod = a_elem * b_elem;
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

	    tmp1 = sum * alpha_i;
	    c_elem = c_i[cij];
	    tmp2 = c_elem * beta_i;
	    tmp1 = tmp1 + tmp2;
	    c_i[cij] = tmp1;
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
      const double *a_i = a;
      const float *b_i = b;

      /* Output Matrix */
      double *c_i = c;

      /* Input Scalars */
      double alpha_i = alpha;
      double beta_i = beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
      float b_elem;
      double c_elem;
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;

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

      FPU_FIX_START;

      /* Ajustment to increments */







      /* alpha = 0.  In this case, just return beta * C */
      if (alpha_i == 0.0) {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci) {
	  cij = ci;
	  for (j = 0; j < n; j++, cij += inccij) {
	    c_elem = c_i[cij];
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = c_elem * split;
	      a1 = con - c_elem;
	      a1 = con - a1;
	      a2 = c_elem - a1;
	      con = beta_i * split;
	      b1 = con - beta_i;
	      b1 = con - b1;
	      b2 = beta_i - b1;

	      head_tmp1 = c_elem * beta_i;
	      tail_tmp1 =
		(((a1 * b1 - head_tmp1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    c_i[cij] = head_tmp1;
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

	      head_sum = tail_sum = 0.0;

	      for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
		a_elem = a_i[aih];
		b_elem = b_i[bhj];
		if (transa == blas_conj_trans) {

		}
		if (transb == blas_conj_trans) {

		}
		{
		  double dt = (double) b_elem;
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = a_elem * split;
		    a1 = con - a_elem;
		    a1 = con - a1;
		    a2 = a_elem - a1;
		    con = dt * split;
		    b1 = con - dt;
		    b1 = con - b1;
		    b2 = dt - b1;

		    head_prod = a_elem * dt;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
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
	      c_i[cij] = head_sum;
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
		{
		  double dt = (double) b_elem;
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = a_elem * split;
		    a1 = con - a_elem;
		    a1 = con - a1;
		    a2 = a_elem - a1;
		    con = dt * split;
		    b1 = con - dt;
		    b1 = con - b1;
		    b2 = dt - b1;

		    head_prod = a_elem * dt;
		    tail_prod =
		      (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		}
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

	      c_elem = c_i[cij];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem * split;
		a1 = con - c_elem;
		a1 = con - a1;
		a2 = c_elem - a1;
		con = beta_i * split;
		b1 = con - beta_i;
		b1 = con - b1;
		b2 = beta_i - b1;

		head_tmp2 = c_elem * beta_i;
		tail_tmp2 =
		  (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_tmp1 = head_sum;
	      tail_tmp1 = tail_sum;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_tmp2 + head_tmp1;
		bv = s1 - head_tmp2;
		s2 = ((head_tmp1 - bv) + (head_tmp2 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_tmp2 + tail_tmp1;
		bv = t1 - tail_tmp2;
		t2 = ((tail_tmp1 - bv) + (tail_tmp2 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_tmp1 = t1 + t2;
		tail_tmp1 = t2 - (head_tmp1 - t1);
	      }
	      c_i[cij] = head_tmp1;
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
	      {
		double dt = (double) b_elem;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = a_elem * split;
		  a1 = con - a_elem;
		  a1 = con - a1;
		  a2 = a_elem - a1;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  head_prod = a_elem * dt;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
	      }
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
	      /* Compute double-double = double-double * double. */
	      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	      con = head_sum * split;
	      a11 = con - head_sum;
	      a11 = con - a11;
	      a21 = head_sum - a11;
	      con = alpha_i * split;
	      b1 = con - alpha_i;
	      b1 = con - b1;
	      b2 = alpha_i - b1;

	      c11 = head_sum * alpha_i;
	      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	      c2 = tail_sum * alpha_i;
	      t1 = c11 + c2;
	      t2 = (c2 - (t1 - c11)) + c21;

	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    c_elem = c_i[cij];
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = c_elem * split;
	      a1 = con - c_elem;
	      a1 = con - a1;
	      a2 = c_elem - a1;
	      con = beta_i * split;
	      b1 = con - beta_i;
	      b1 = con - b1;
	      b2 = beta_i - b1;

	      head_tmp2 = c_elem * beta_i;
	      tail_tmp2 =
		(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_tmp1 + head_tmp2;
	      bv = s1 - head_tmp1;
	      s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_tmp1 + tail_tmp2;
	      bv = t1 - tail_tmp1;
	      t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    c_i[cij] = head_tmp1;
	  }
	}

      }

      FPU_FIX_STOP;

      break;
    }
  }
}
