#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void BLAS_ssymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, float *alpha,
			int alpha_flag, float *beta, int beta_flag, float *a,
			int lda, float *b, int ldb, float *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_ssymm{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) float*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) float*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) float*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) float*
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  float c_elem;
  float a_elem;
  float b_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *b_vec;

  float *c_i = c;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *a_i = a;
  float *b_i = b;

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





  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_sdot_testgen(m_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, &c_elem,
		      &head_r_true_elem, &tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	sge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	sge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_sdot_testgen(m_i, i, m_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			&c_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem;
      head_r_true[cij] = head_r_true_elem;
      tail_r_true[cij] = tail_r_true_elem;
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem = xrand(seed);
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = xrand(seed);
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_sdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &head_r_true_elem, &tail_r_true_elem);

	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true_elem;
	tail_r_true[cij] = tail_r_true_elem;
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dsymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, double *alpha,
			int alpha_flag, double *beta, int beta_flag,
			double *a, int lda, double *b, int ldb, double *c,
			int ldc, int *seed, double *head_r_true,
			double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymm{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) double*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) double*
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double a_elem;
  double b_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  double *b_vec;

  double *c_i = c;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  double *b_i = b;

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





  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_ddot_testgen(m_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, &c_elem,
		      &head_r_true_elem, &tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	dge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	dge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_ddot_testgen(m_i, i, m_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			&c_elem, &head_r_true_elem, &tail_r_true_elem);

      dsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem;
      head_r_true[cij] = head_r_true_elem;
      tail_r_true[cij] = tail_r_true_elem;
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem = xrand(seed);
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = xrand(seed);
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_ddot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &head_r_true_elem, &tail_r_true_elem);

	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true_elem;
	tail_r_true[cij] = tail_r_true_elem;
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_csymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymm{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  float c_elem[2];
  float a_elem[2];
  float b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *b_vec;

  float *c_i = (float *) c;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_cdot_testgen(m_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, c_elem,
		      head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	cge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	cge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_cdot_testgen(m_i, i, m_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			c_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = xrand(seed);
	b_elem[1] = xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
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
void BLAS_zsymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymm{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  double a_elem[2];
  double b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *b_i = (double *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_testgen(m_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, c_elem,
		      head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	zge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	zge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_testgen(m_i, i, m_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			c_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = xrand(seed);
	b_elem[1] = xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  zge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  zge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
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
void BLAS_csymm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, float *a, int lda,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymm_s_s{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 * a       (input/output) float*
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  float c_elem[2];
  float a_elem;
  float b_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

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

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_cdot_s_s_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	sge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	sge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_cdot_s_s_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      ssy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {

    float *aa_vec;
    float *bb_vec;

    aa_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
    if (m_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    bb_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = (float) xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      {
	int r;
	for (r = 0; r < m_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	{
	  int r;
	  for (r = 0; r < m_i; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}


	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  aa_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(aa_vec);
    blas_free(bb_vec);
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_csymm_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, float *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymm_s_c{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 * a       (input/output) float*
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  float c_elem[2];
  float a_elem;
  float b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *b_vec;

  float *c_i = (float *) c;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *b_i = (float *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;

  incb *= 2;
  a_vec = (float *) blas_malloc(m_i * sizeof(float));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_cdot_c_s_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	cge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	cge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_cdot_c_s_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      ssy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {

    float *aa_vec;


    aa_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
    if (m_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      {
	int r;
	for (r = 0; r < m_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  b_vec,
			  aa_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(aa_vec);

  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_csymm_c_s_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_csymm_c_s{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  float c_elem[2];
  float a_elem[2];
  float b_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *b_vec;

  float *c_i = (float *) c;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_cdot_s_c_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	sge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	sge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_cdot_s_c_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {


    float *bb_vec;


    bb_vec = (float *) blas_malloc(m_i * sizeof(float) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }



    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	{
	  int r;
	  for (r = 0; r < m_i; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}


	BLAS_cdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  a_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

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
void BLAS_zsymm_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, double *a, int lda,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymm_d_d{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 * a       (input/output) double*
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  double a_elem;
  double b_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

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

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_d_d_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	dge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	dge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_d_d_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      dsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {

    double *aa_vec;
    double *bb_vec;

    aa_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
    if (m_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    bb_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = (float) xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      {
	int r;
	for (r = 0; r < m_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	{
	  int r;
	  for (r = 0; r < m_i; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}


	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  aa_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(aa_vec);
    blas_free(bb_vec);
  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zsymm_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, double *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymm_d_z{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 * a       (input/output) double*
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  double a_elem;
  double b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *b_i = (double *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  inca = incb = 1;

  incb *= 2;
  a_vec = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < m_i * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_z_d_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	zge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	zge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_z_d_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      dsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {

    double *aa_vec;


    aa_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
    if (m_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      {
	int r;
	for (r = 0; r < m_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  zge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  zge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  b_vec,
			  aa_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true_elem[0];
	head_r_true[cij + 1] = head_r_true_elem[1];
	tail_r_true[cij] = tail_r_true_elem[0];
	tail_r_true[cij + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(aa_vec);

  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zsymm_z_d_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zsymm_z_d{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  double a_elem[2];
  double b_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *b_i = b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_d_z_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	dge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	dge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_d_z_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {


    double *bb_vec;


    bb_vec = (double *) blas_malloc(m_i * sizeof(double) * 2);
    if (m_i > 0 && bb_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    if (alpha_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = xrand(seed);
      c_elem[1] = xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }



    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	{
	  int r;
	  for (r = 0; r < m_i; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}


	BLAS_zdot_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  a_vec,
			  seed, c_elem, head_r_true_elem, tail_r_true_elem);

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
void BLAS_dsymm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, float *a, int lda,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymm_s_s{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) float*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) double*
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  float a_elem;
  float b_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *b_vec;

  double *c_i = c;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  float *b_i = b;

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





  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_ddot_s_s_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	sge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	sge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_ddot_s_s_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem;
      head_r_true[cij] = head_r_true_elem;
      tail_r_true[cij] = tail_r_true_elem;
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem = (float) xrand(seed);
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) xrand(seed);
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = (float) xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_ddot_s_s_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &head_r_true_elem, &tail_r_true_elem);

	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true_elem;
	tail_r_true[cij] = tail_r_true_elem;
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dsymm_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, float *a, int lda,
			    double *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymm_s_d{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) double*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) double*
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  float a_elem;
  double b_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  double *b_vec;

  double *c_i = c;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  double *b_i = b;

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





  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_ddot_d_s_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	dge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	dge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_ddot_d_s_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem;
      head_r_true[cij] = head_r_true_elem;
      tail_r_true[cij] = tail_r_true_elem;
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem = (float) xrand(seed);
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) xrand(seed);
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = (float) xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      ssy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  dge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  dge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_ddot_d_s_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &head_r_true_elem, &tail_r_true_elem);

	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true_elem;
	tail_r_true[cij] = tail_r_true_elem;
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_dsymm_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, double *a, int lda,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymm_d_s{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) float*
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) double*
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem;
  double a_elem;
  float b_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  float *b_vec;

  double *c_i = c;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  float *b_i = b;

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





  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_ddot_s_d_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem;
    head_r_true[cij] = head_r_true_elem;
    tail_r_true[cij] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	sge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	sge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_ddot_s_d_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &head_r_true_elem, &tail_r_true_elem);

      dsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem;
      head_r_true[cij] = head_r_true_elem;
      tail_r_true[cij] = tail_r_true_elem;
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem = (float) xrand(seed);
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) xrand(seed);
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }






    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem = (float) xrand(seed);
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      dsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  sge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  sge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_ddot_s_d_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &head_r_true_elem, &tail_r_true_elem);

	c_i[cij] = c_elem;
	head_r_true[cij] = head_r_true_elem;
	tail_r_true[cij] = tail_r_true_elem;
      }
    }



  }

  blas_free(a_vec);
  blas_free(b_vec);
}
void BLAS_zsymm_c_c_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zsymm_c_c{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  float a_elem[2];
  float b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_c_c_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	cge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	cge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_c_c_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_zdot_c_c_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
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
void BLAS_zsymm_c_z_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zsymm_c_z{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  float a_elem[2];
  double b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  double *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  double *b_i = (double *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_z_c_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	zge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	zge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_z_c_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      csy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  zge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  zge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_zdot_z_c_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
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
void BLAS_zsymm_z_c_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zsymm_z_c{_x}
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
 *           which half of the symmetric matrix a is to be stored.
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
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
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
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *head_r_true
 *         the leading part of the truth in double-double.
 *
 * double  (output) *tail_r_true
 *         the trailing part of the truth in double-double
 *
 */
{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  double c_elem[2];
  double a_elem[2];
  float b_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  float *b_vec;

  double *c_i = (double *) c;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  float *b_i = (float *) b;

  if (side == blas_left_side) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

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

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */

    BLAS_zdot_c_z_testgen(m_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  head_r_true_elem, tail_r_true_elem);

    cij = 0;
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    head_r_true[cij] = head_r_true_elem[0];
    head_r_true[cij + 1] = head_r_true_elem[1];
    tail_r_true[cij] = tail_r_true_elem[0];
    tail_r_true[cij + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
	cge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
	cge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      BLAS_zdot_c_z_testgen(m_i, i, m_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      head_r_true[cij] = head_r_true_elem[0];
      head_r_true[cij + 1] = head_r_true_elem[1];
      tail_r_true[cij] = tail_r_true_elem[0];
      tail_r_true[cij + 1] = tail_r_true_elem[1];
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	head_r_true[cij] = head_r_true[ci];
	tail_r_true[cij] = tail_r_true[ci];
	head_r_true[cij + 1] = head_r_true[ci + 1];
	tail_r_true[cij + 1] = tail_r_true[ci + 1];
      }
    }
  } else {







    if (alpha_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      c_elem[0] = (float) xrand(seed);
      c_elem[1] = (float) xrand(seed);
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
	(order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
	b_elem[0] = (float) xrand(seed);
	b_elem[1] = (float) xrand(seed);
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      zsy_copy_row(order, uplo, m_i, a, lda, a_vec, i);


      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

	if (side == blas_left_side)
	  cge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
	else
	  cge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);



	BLAS_zdot_c_z_testgen(m_i, m_i, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
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
