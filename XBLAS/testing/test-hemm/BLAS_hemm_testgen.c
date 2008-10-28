#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void BLAS_chemm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_chemm{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) void*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  float *a1;
  float *a2;
  float *c1;
  float *c2;
  float *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  float *a_vec;
  float *b_vec;

  float *c_i = (float *) c;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;
  a_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    a1 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_ssymm_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_sskew_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a complex matrix.  Since B is generated
       as a real matrix, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    B       beta    C    R (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale R by 1+i, since they are
       computed in double-double precision.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by 1+i. */
      alpha_i[1] = 0.0;
      beta_i[1] = beta_i[0];
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by 1+i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
    } else {
      ab = 3;
      /* multiply alpha by 1+i, beta by 2i. */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = 2.0 * beta_i[0];
      beta_i[0] = 0.0;
    }


    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	if (ab == 0) {
	  b_i[bij] = 0.0;
	  b_i[bij + 1] = b0[mij];
	} else {
	  b_i[bij] = b0[mij];
	  b_i[bij + 1] = b0[mij];
	}
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 0) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else if (ab == 2) {
	  c_i[cij] = -2.0 * c2[mij];
	  c_i[cij + 1] = 2.0 * c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	head_r_elem1 = head_r1_true[mij];
	tail_r_elem1 = tail_r1_true[mij];

	head_r_elem2 = head_r2_true[mij];
	tail_r_elem2 = tail_r2_true[mij];

	if (ab == 0) {
	  head_r_true[cij] = -head_r_elem2;
	  tail_r_true[cij] = -tail_r_elem2;
	  head_r_true[cij + 1] = head_r_elem1;
	  tail_r_true[cij + 1] = tail_r_elem1;
	} else if (ab == 1) {
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	} else {

	  /* Real part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem2 * split;
	    a11 = con - head_r_elem2;
	    a11 = con - a11;
	    a21 = head_r_elem2 - a11;
	    con = -2.0 * split;
	    b1 = con - -2.0;
	    b1 = con - b1;
	    b2 = -2.0 - b1;

	    c11 = head_r_elem2 * -2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem2 * -2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij] = head_r_elem;
	  tail_r_true[cij] = tail_r_elem;

	  /* Imaginary Part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem1 * split;
	    a11 = con - head_r_elem1;
	    a11 = con - a11;
	    a21 = head_r_elem1 - a11;
	    con = 2.0 * split;
	    b1 = con - 2.0;
	    b1 = con - b1;
	    b2 = 2.0 - b1;

	    c11 = head_r_elem1 * 2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem1 * 2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij + 1] = head_r_elem;
	  tail_r_true[cij + 1] = tail_r_elem;
	}
      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    float a_elem[2];
    float b_elem[2];
    float c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */


    if (alpha_flag == 0) {
      alpha_i[0] = xrand(seed);
      alpha_i[1] = xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = xrand(seed);
      beta_i[1] = xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = xrand(seed);
	b_elem[1] = xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      che_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec, b_vec, seed,
			  c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zhemm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemm{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) void*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  double *a1;
  double *a2;
  double *c1;
  double *c2;
  double *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  double *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;
  a_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    a1 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_dsymm_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_dskew_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a complex matrix.  Since B is generated
       as a real matrix, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    B       beta    C    R (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale R by 1+i, since they are
       computed in double-double precision.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by 1+i. */
      alpha_i[1] = 0.0;
      beta_i[1] = beta_i[0];
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by 1+i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
    } else {
      ab = 3;
      /* multiply alpha by 1+i, beta by 2i. */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = 2.0 * beta_i[0];
      beta_i[0] = 0.0;
    }


    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	if (ab == 0) {
	  b_i[bij] = 0.0;
	  b_i[bij + 1] = b0[mij];
	} else {
	  b_i[bij] = b0[mij];
	  b_i[bij + 1] = b0[mij];
	}
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 0) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else if (ab == 2) {
	  c_i[cij] = -2.0 * c2[mij];
	  c_i[cij + 1] = 2.0 * c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	head_r_elem1 = head_r1_true[mij];
	tail_r_elem1 = tail_r1_true[mij];

	head_r_elem2 = head_r2_true[mij];
	tail_r_elem2 = tail_r2_true[mij];

	if (ab == 0) {
	  head_r_true[cij] = -head_r_elem2;
	  tail_r_true[cij] = -tail_r_elem2;
	  head_r_true[cij + 1] = head_r_elem1;
	  tail_r_true[cij + 1] = tail_r_elem1;
	} else if (ab == 1) {
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	} else {

	  /* Real part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem2 * split;
	    a11 = con - head_r_elem2;
	    a11 = con - a11;
	    a21 = head_r_elem2 - a11;
	    con = -2.0 * split;
	    b1 = con - -2.0;
	    b1 = con - b1;
	    b2 = -2.0 - b1;

	    c11 = head_r_elem2 * -2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem2 * -2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij] = head_r_elem;
	  tail_r_true[cij] = tail_r_elem;

	  /* Imaginary Part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem1 * split;
	    a11 = con - head_r_elem1;
	    a11 = con - a11;
	    a21 = head_r_elem1 - a11;
	    con = 2.0 * split;
	    b1 = con - 2.0;
	    b1 = con - b1;
	    b2 = 2.0 - b1;

	    c11 = head_r_elem1 * 2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem1 * 2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij + 1] = head_r_elem;
	  tail_r_true[cij + 1] = tail_r_elem;
	}
      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    double a_elem[2];
    double b_elem[2];
    double c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */


    if (alpha_flag == 0) {
      alpha_i[0] = xrand(seed);
      alpha_i[1] = xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = xrand(seed);
      beta_i[1] = xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = xrand(seed);
	b_elem[1] = xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zhe_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  zge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  zge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec, b_vec, seed,
			  c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zhemm_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemm_c_z{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) void*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  float *a1;
  float *a2;
  double *c1;
  double *c2;
  double *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  float *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = a;
  double *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;
  a_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    a1 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_dsymm_s_d_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_dskew_s_d_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a complex matrix.  Since B is generated
       as a real matrix, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    B       beta    C    R (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale R by 1+i, since they are
       computed in double-double precision.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by 1+i. */
      alpha_i[1] = 0.0;
      beta_i[1] = beta_i[0];
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by 1+i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
    } else {
      ab = 3;
      /* multiply alpha by 1+i, beta by 2i. */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = 2.0 * beta_i[0];
      beta_i[0] = 0.0;
    }


    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	if (ab == 0) {
	  b_i[bij] = 0.0;
	  b_i[bij + 1] = b0[mij];
	} else {
	  b_i[bij] = b0[mij];
	  b_i[bij + 1] = b0[mij];
	}
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 0) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else if (ab == 2) {
	  c_i[cij] = -2.0 * c2[mij];
	  c_i[cij + 1] = 2.0 * c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	head_r_elem1 = head_r1_true[mij];
	tail_r_elem1 = tail_r1_true[mij];

	head_r_elem2 = head_r2_true[mij];
	tail_r_elem2 = tail_r2_true[mij];

	if (ab == 0) {
	  head_r_true[cij] = -head_r_elem2;
	  tail_r_true[cij] = -tail_r_elem2;
	  head_r_true[cij + 1] = head_r_elem1;
	  tail_r_true[cij + 1] = tail_r_elem1;
	} else if (ab == 1) {
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	} else {

	  /* Real part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem2 * split;
	    a11 = con - head_r_elem2;
	    a11 = con - a11;
	    a21 = head_r_elem2 - a11;
	    con = -2.0 * split;
	    b1 = con - -2.0;
	    b1 = con - b1;
	    b2 = -2.0 - b1;

	    c11 = head_r_elem2 * -2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem2 * -2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij] = head_r_elem;
	  tail_r_true[cij] = tail_r_elem;

	  /* Imaginary Part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem1 * split;
	    a11 = con - head_r_elem1;
	    a11 = con - a11;
	    a21 = head_r_elem1 - a11;
	    con = 2.0 * split;
	    b1 = con - 2.0;
	    b1 = con - b1;
	    b2 = 2.0 - b1;

	    c11 = head_r_elem1 * 2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem1 * 2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij + 1] = head_r_elem;
	  tail_r_true[cij + 1] = tail_r_elem;
	}
      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    float a_elem[2];
    double b_elem[2];
    double c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */


    if (alpha_flag == 0) {
      alpha_i[0] = (float) xrand(seed);
      alpha_i[1] = (float) xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = (float) xrand(seed);
      beta_i[1] = (float) xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      che_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  zge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  zge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	BLAS_zdot_c_z_testgen(m_i, m_i, 0, norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec, b_vec, seed,
			      c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zhemm_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemm_z_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) void*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  double *a1;
  double *a2;
  double *c1;
  double *c2;
  float *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  double *a_vec;
  float *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  float *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;
  a_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    a1 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_dsymm_d_s_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_dskew_d_s_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a complex matrix.  Since B is generated
       as a real matrix, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    B       beta    C    R (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale R by 1+i, since they are
       computed in double-double precision.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by 1+i. */
      alpha_i[1] = 0.0;
      beta_i[1] = beta_i[0];
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by 1+i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
    } else {
      ab = 3;
      /* multiply alpha by 1+i, beta by 2i. */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = 2.0 * beta_i[0];
      beta_i[0] = 0.0;
    }


    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	if (ab == 0) {
	  b_i[bij] = 0.0;
	  b_i[bij + 1] = b0[mij];
	} else {
	  b_i[bij] = b0[mij];
	  b_i[bij + 1] = b0[mij];
	}
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 0) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else if (ab == 2) {
	  c_i[cij] = -2.0 * c2[mij];
	  c_i[cij + 1] = 2.0 * c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	head_r_elem1 = head_r1_true[mij];
	tail_r_elem1 = tail_r1_true[mij];

	head_r_elem2 = head_r2_true[mij];
	tail_r_elem2 = tail_r2_true[mij];

	if (ab == 0) {
	  head_r_true[cij] = -head_r_elem2;
	  tail_r_true[cij] = -tail_r_elem2;
	  head_r_true[cij + 1] = head_r_elem1;
	  tail_r_true[cij + 1] = tail_r_elem1;
	} else if (ab == 1) {
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	} else {

	  /* Real part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem2 * split;
	    a11 = con - head_r_elem2;
	    a11 = con - a11;
	    a21 = head_r_elem2 - a11;
	    con = -2.0 * split;
	    b1 = con - -2.0;
	    b1 = con - b1;
	    b2 = -2.0 - b1;

	    c11 = head_r_elem2 * -2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem2 * -2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij] = head_r_elem;
	  tail_r_true[cij] = tail_r_elem;

	  /* Imaginary Part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem1 * split;
	    a11 = con - head_r_elem1;
	    a11 = con - a11;
	    a21 = head_r_elem1 - a11;
	    con = 2.0 * split;
	    b1 = con - 2.0;
	    b1 = con - b1;
	    b2 = 2.0 - b1;

	    c11 = head_r_elem1 * 2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem1 * 2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij + 1] = head_r_elem;
	  tail_r_true[cij + 1] = tail_r_elem;
	}
      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    double a_elem[2];
    float b_elem[2];
    double c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */


    if (alpha_flag == 0) {
      alpha_i[0] = (float) xrand(seed);
      alpha_i[1] = (float) xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = (float) xrand(seed);
      beta_i[1] = (float) xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zhe_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	BLAS_zdot_z_c_testgen(m_i, m_i, 0, norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec, b_vec, seed,
			      c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zhemm_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemm_c_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) void*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  float *a1;
  float *a2;
  double *c1;
  double *c2;
  float *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  float *a_vec;
  float *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = a;
  float *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;
  a_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    a1 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_dsymm_s_s_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_dskew_s_s_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a complex matrix.  Since B is generated
       as a real matrix, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    B       beta    C    R (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale R by 1+i, since they are
       computed in double-double precision.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by 1+i. */
      alpha_i[1] = 0.0;
      beta_i[1] = beta_i[0];
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by 1+i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
    } else {
      ab = 3;
      /* multiply alpha by 1+i, beta by 2i. */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = 2.0 * beta_i[0];
      beta_i[0] = 0.0;
    }


    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	if (ab == 0) {
	  b_i[bij] = 0.0;
	  b_i[bij + 1] = b0[mij];
	} else {
	  b_i[bij] = b0[mij];
	  b_i[bij + 1] = b0[mij];
	}
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 0) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else if (ab == 2) {
	  c_i[cij] = -2.0 * c2[mij];
	  c_i[cij + 1] = 2.0 * c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	head_r_elem1 = head_r1_true[mij];
	tail_r_elem1 = tail_r1_true[mij];

	head_r_elem2 = head_r2_true[mij];
	tail_r_elem2 = tail_r2_true[mij];

	if (ab == 0) {
	  head_r_true[cij] = -head_r_elem2;
	  tail_r_true[cij] = -tail_r_elem2;
	  head_r_true[cij + 1] = head_r_elem1;
	  tail_r_true[cij + 1] = tail_r_elem1;
	} else if (ab == 1) {
	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	} else {

	  /* Real part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem2 * split;
	    a11 = con - head_r_elem2;
	    a11 = con - a11;
	    a21 = head_r_elem2 - a11;
	    con = -2.0 * split;
	    b1 = con - -2.0;
	    b1 = con - b1;
	    b2 = -2.0 - b1;

	    c11 = head_r_elem2 * -2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem2 * -2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij] = head_r_elem;
	  tail_r_true[cij] = tail_r_elem;

	  /* Imaginary Part */
	  {
	    /* Compute double-double = double-double * double. */
	    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = head_r_elem1 * split;
	    a11 = con - head_r_elem1;
	    a11 = con - a11;
	    a21 = head_r_elem1 - a11;
	    con = 2.0 * split;
	    b1 = con - 2.0;
	    b1 = con - b1;
	    b2 = 2.0 - b1;

	    c11 = head_r_elem1 * 2.0;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	    c2 = tail_r_elem1 * 2.0;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;

	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }
	  head_r_true[cij + 1] = head_r_elem;
	  tail_r_true[cij + 1] = tail_r_elem;
	}
      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    float a_elem[2];
    float b_elem[2];
    double c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */


    if (alpha_flag == 0) {
      alpha_i[0] = (float) xrand(seed);
      alpha_i[1] = (float) xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = (float) xrand(seed);
      beta_i[1] = (float) xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      che_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	BLAS_zdot_c_c_testgen(m_i, m_i, 0, norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec, b_vec, seed,
			      c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}

void BLAS_zhemm_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemm_z_d{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) double*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  double *a1;
  double *a2;
  double *c1;
  double *c2;
  double *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  double *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;



  inca = incb = 1;
  inca *= 2;

  a_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if (randomize == 0) {
    a1 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(m_i * m_i * sizeof(double));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_dsymm_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_dskew_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a real matrix. 

       There are four cases to consider, depending on the 
       values of alpha and beta.

       values                             scaling
       alpha  beta         alpha    A    B    beta    C     R (truth)
       0    1      1            
       1    1      ?                              -i     i     
       2    ?      1            i                        i     i
       3    ?      ?           1+i                1+i         1+i

       Note that we can afford to scale truth by (1+i) since they
       are computed in double-double.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by -i. */
      alpha_i[1] = 0.0;
      beta_i[1] = -beta_i[0];
      beta_i[0] = 0.0;
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
      alpha_i[0] = 0.0;
    } else {
      ab = 3;
      /* multiply alpha, beta by (1 + i). */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = beta_i[0];
    }

    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	b_i[bij] = b0[mij];
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 1 || ab == 2) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	if (ab == 0 || ab == 1) {
	  head_r_true[cij] = head_r1_true[mij];
	  tail_r_true[cij] = tail_r1_true[mij];
	  head_r_true[cij + 1] = head_r2_true[mij];
	  tail_r_true[cij + 1] = tail_r2_true[mij];
	} else if (ab == 2) {
	  head_r_true[cij] = -head_r2_true[mij];
	  tail_r_true[cij] = -tail_r2_true[mij];
	  head_r_true[cij + 1] = head_r1_true[mij];
	  tail_r_true[cij + 1] = tail_r1_true[mij];
	} else {
	  head_r_elem1 = head_r1_true[mij];
	  tail_r_elem1 = tail_r1_true[mij];

	  head_r_elem2 = head_r2_true[mij];
	  tail_r_elem2 = tail_r2_true[mij];

	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	}

      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    double a_elem[2];
    double b_elem;
    double c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */
    double *bb_vec;
    bb_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      alpha_i[0] = xrand(seed);
      alpha_i[1] = xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = xrand(seed);
      beta_i[1] = xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zhe_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	{
	  int k;
	  for (k = 0; k < m_i; k++) {
	    bb_vec[2 * k] = b_vec[k];
	    bb_vec[2 * k + 1] = 0.0;
	  }
	}
	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec, bb_vec, seed,
			  c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(bb_vec);

  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_chemm_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_chemm_c_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
 * 
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) float*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) void*
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * head_r_true  (output) double *
 *         the leading part of the truth in double-double.
 *
 * tail_r_true  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{

  /* Strategy:  
     R1 = alpha * A1 * B + beta * C1
     R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
   */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  float *a1;
  float *a2;
  float *c1;
  float *c2;
  float *b0;

  double *head_r1_true, *tail_r1_true;
  double *head_r2_true, *tail_r2_true;

  double head_r_elem1, tail_r_elem1;
  double head_r_elem2, tail_r_elem2;
  double head_r_elem, tail_r_elem;

  float *a_vec;
  float *b_vec;

  float *c_i = (float *) c;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  incci *= 2;
  inccij *= 2;
  incai *= 2;
  incaij *= 2;



  inca = incb = 1;
  inca *= 2;

  a_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if (randomize == 0) {
    a1 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(m_i * m_i * sizeof(float));
    if (m_i * m_i > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * m_i; ++i) {
      a1[i] = a2[i] = 0.0;
    }
    c1 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && c1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    c2 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && c2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    b0 = (float *) blas_malloc(m_i * n_i * sizeof(float));
    if (m_i * n_i > 0 && b0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < m_i * n_i; ++i) {
      c1[i] = c2[i] = b0[i] = 0.0;
    }
    head_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r1_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    tail_r2_true = (double *) blas_malloc(m_i * n_i * sizeof(double));
    if (m_i * n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of matrix A, and matrix B.
       Note that Re(A) is a symmetric matrix. */
    BLAS_ssymm_testgen
      (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i,
       beta_flag, a1, m_i, b0, ld, c1, ld, seed, head_r1_true, tail_r1_true);

    BLAS_sskew_testgen
      (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld,
       c2, ld, seed, head_r2_true, tail_r2_true);


    /* The case where B is a real matrix. 

       There are four cases to consider, depending on the 
       values of alpha and beta.

       values                             scaling
       alpha  beta         alpha    A    B    beta    C     R (truth)
       0    1      1            
       1    1      ?                              -i     i     
       2    ?      1            i                        i     i
       3    ?      ?           1+i                1+i         1+i

       Note that we can afford to scale truth by (1+i) since they
       are computed in double-double.
     */

    if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
      ab = 0;
      alpha_i[1] = beta_i[1] = 0.0;	/* set alpha, beta to be 1. */
    } else if (alpha_i[0] == 1.0) {
      ab = 1;
      /* set alpha to 1, multiply beta by -i. */
      alpha_i[1] = 0.0;
      beta_i[1] = -beta_i[0];
      beta_i[0] = 0.0;
    } else if (beta_i[0] == 1.0) {
      ab = 2;
      /* set beta to 1, multiply alpha by i. */
      beta_i[1] = 0.0;
      alpha_i[1] = alpha_i[0];
      alpha_i[0] = 0.0;
    } else {
      ab = 3;
      /* multiply alpha, beta by (1 + i). */
      alpha_i[1] = alpha_i[0];
      beta_i[1] = beta_i[0];
    }

    /* Now fill in a */
    for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < m_i;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in b */
    for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
      for (j = 0, bij = bi, mij = mi; j < n_i;
	   j++, bij += incbij, mij += incij) {
	b_i[bij] = b0[mij];
      }
    }

    /* Fill in c */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {
	if (ab == 1 || ab == 2) {
	  c_i[cij] = -c2[mij];
	  c_i[cij + 1] = c1[mij];
	} else {
	  c_i[cij] = c1[mij];
	  c_i[cij + 1] = c2[mij];
	}
      }
    }

    /* Fill in truth */
    for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
      for (j = 0, cij = ci, mij = mi; j < n_i;
	   j++, cij += inccij, mij += incij) {

	if (ab == 0 || ab == 1) {
	  head_r_true[cij] = head_r1_true[mij];
	  tail_r_true[cij] = tail_r1_true[mij];
	  head_r_true[cij + 1] = head_r2_true[mij];
	  tail_r_true[cij + 1] = tail_r2_true[mij];
	} else if (ab == 2) {
	  head_r_true[cij] = -head_r2_true[mij];
	  tail_r_true[cij] = -tail_r2_true[mij];
	  head_r_true[cij + 1] = head_r1_true[mij];
	  tail_r_true[cij + 1] = tail_r1_true[mij];
	} else {
	  head_r_elem1 = head_r1_true[mij];
	  tail_r_elem1 = tail_r1_true[mij];

	  head_r_elem2 = head_r2_true[mij];
	  tail_r_elem2 = tail_r2_true[mij];

	  {
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_r_elem1 + head_r_elem2;
	    bv = s1 - head_r_elem1;
	    s2 = ((head_r_elem2 - bv) + (head_r_elem1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_r_elem1 + tail_r_elem2;
	    bv = t1 - tail_r_elem1;
	    t2 = ((tail_r_elem2 - bv) + (tail_r_elem1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_r_elem = t1 + t2;
	    tail_r_elem = t2 - (head_r_elem - t1);
	  }

	  /* Set the imaginary part to  R1 + R2 */
	  tail_r_true[cij + 1] = tail_r_elem;
	  head_r_true[cij + 1] = head_r_elem;

	  /* Set the real part to R1 - R2. */
	  {
	    double head_bt, tail_bt;
	    head_bt = -head_r_elem2;
	    tail_bt = -tail_r_elem2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_r_elem1 + head_bt;
	      bv = s1 - head_r_elem1;
	      s2 = ((head_bt - bv) + (head_r_elem1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_r_elem1 + tail_bt;
	      bv = t1 - tail_r_elem1;
	      t2 = ((tail_bt - bv) + (tail_r_elem1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_r_elem = t1 + t2;
	      tail_r_elem = t2 - (head_r_elem - t1);
	    }
	  }
	  tail_r_true[cij] = tail_r_elem;
	  head_r_true[cij] = head_r_elem;
	}

      }
    }

    blas_free(a1);
    blas_free(a2);
    blas_free(c1);
    blas_free(c2);
    blas_free(b0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);

  } else {
    /* random matrix */
    float a_elem[2];
    float b_elem;
    float c_elem[2];
    double head_r_true_elem[2], tail_r_true_elem[2];

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */
    float *bb_vec;
    bb_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      alpha_i[0] = xrand(seed);
      alpha_i[1] = xrand(seed);
    }
    if (beta_flag == 0) {
      beta_i[0] = xrand(seed);
      beta_i[1] = xrand(seed);
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
	if (i == j)
	  a_i[aij + 1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      che_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);

	/* copy the real b_vec into complex bb_vec, so that 
	   pure complex test case generator can be called. */

	{
	  int k;
	  for (k = 0; k < m_i; k++) {
	    bb_vec[2 * k] = b_vec[k];
	    bb_vec[2 * k + 1] = 0.0;
	  }
	}
	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec, bb_vec, seed,
			  c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(bb_vec);

  }

  blas_free(a_vec);
  blas_free(b_vec);
}

void BLAS_sskew_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, float *alpha, int alpha_flag,
			float *beta, int beta_flag, float *a, int lda,
			float *b, int ldb, float *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  float c_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *b_vec;

  float *c_i = c;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;


  a_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (side == blas_left_side)
    sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0);

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    sskew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    BLAS_sdot_testgen(m_i, i, m_i - i, norm,
		      blas_no_conj, alpha, 1,
		      beta, 1, b_vec, a_vec, seed,
		      &c_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;
  }

  /* Now fill in c and r_true */

  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      c_i[cij] = c_i[ci];
      head_r_true[cij] = head_r_true[ci];
      tail_r_true[cij] = tail_r_true[ci];
    }
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dskew_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, double *alpha, int alpha_flag,
			double *beta, int beta_flag, double *a, int lda,
			double *b, int ldb, double *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  double *b_vec;

  double *c_i = c;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;


  a_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (side == blas_left_side)
    dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0);

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    dskew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    BLAS_ddot_testgen(m_i, i, m_i - i, norm,
		      blas_no_conj, alpha, 1,
		      beta, 1, b_vec, a_vec, seed,
		      &c_elem, &head_r_true_elem, &tail_r_true_elem);

    dskew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;
  }

  /* Now fill in c and r_true */

  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      c_i[cij] = c_i[ci];
      head_r_true[cij] = head_r_true[ci];
      tail_r_true[cij] = tail_r_true[ci];
    }
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dskew_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, int lda, float *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  float *b_vec;

  double *c_i = c;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;


  a_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (side == blas_left_side)
    sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0);

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    dskew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    BLAS_ddot_s_d_testgen(m_i, i, m_i - i, norm,
			  blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &head_r_true_elem, &tail_r_true_elem);

    dskew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;
  }

  /* Now fill in c and r_true */

  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      c_i[cij] = c_i[ci];
      head_r_true[cij] = head_r_true[ci];
      tail_r_true[cij] = tail_r_true[ci];
    }
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dskew_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, double *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  double *b_vec;

  double *c_i = c;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;


  a_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (side == blas_left_side)
    dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0);

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    sskew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    BLAS_ddot_d_s_testgen(m_i, i, m_i - i, norm,
			  blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;
  }

  /* Now fill in c and r_true */

  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      c_i[cij] = c_i[ci];
      head_r_true[cij] = head_r_true[ci];
      tail_r_true[cij] = tail_r_true[ci];
    }
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dskew_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, float *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *b_vec;

  double *c_i = c;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;


  a_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (side == blas_left_side)
    sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0);

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    sskew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    BLAS_ddot_s_s_testgen(m_i, i, m_i - i, norm,
			  blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;
  }

  /* Now fill in c and r_true */

  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      c_i[cij] = c_i[ci];
      head_r_true[cij] = head_r_true[ci];
      tail_r_true[cij] = tail_r_true[ci];
    }
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
