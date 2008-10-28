#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

#ifndef ALPHA_USE_IS_ALPHA
#define ALPHA_USE_IS_ALPHA 1
#define ALPHA_USE_IS_BETA 0
#define ALPHA_USE_IS_EITHER -1
#endif

void BLAS_sge_sum_mv_testgen(int norm, enum blas_order_type order,
			     int m, int n, int randomize,
			     float *alpha, int alpha_flag, float *beta,
			     int beta_flag, float *a, int lda, float *b,
			     int ldb, float *x, int incx,
			     float *alpha_use_ptr, float *a_use, float *b_use,
			     int *seed, double *head_r_true,
			     double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_sge_sum_mv{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) float*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) float*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) float*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  float y_elem;
  float beta_zero_fake;
  float a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;
  float multiplier;
  float divider;
  float alpha_use;

  float *a_vec;
  float *x_vec;

  float *alpha_use_ptr_i = alpha_use_ptr;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *a_i = a;
  float *b_i = b;
  float *a_use_i = a_use;
  float *b_use_i = b_use;
  float *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if ((*alpha_i) == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if ((*beta_i) == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if ((*beta_i) == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem = xrand(seed);
      beta_i[0] = y_elem;
    }
    alpha_use = (*beta_i);
  } else {
    if (!alpha_flag) {
      y_elem = xrand(seed);
      alpha_i[0] = y_elem;
    }
    alpha_use = (*alpha_i);
  }
  /* put in return value */
  (*alpha_use_ptr_i) = alpha_use;

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_sdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*alpha_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_sdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*beta_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }




    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem = 0.0;
      BLAS_sdot_testgen(2 * n_i, 0, 2 * n_i, norm,
			blas_no_conj, &alpha_use, 1,
			&beta_zero_fake, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }


  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem = 0.0;

	  BLAS_sdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem = 0.0;

	  BLAS_sdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem = 0.0;
	BLAS_sdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec, a_vec, seed,
			  &y_elem, &head_r_true_elem, &tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem;
	tail_r_true[ri] = tail_r_true_elem;
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  (*alpha_use_ptr_i) = alpha_use;


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {
    {
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	  case 3:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
    case 3:
      (*beta_i) = alpha_use;
      (*alpha_i) = (*beta_i) * multiplier;
      break;
    case 2:			/*should not happen */
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {
      {
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	    case 3:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
      case 3:
	(*alpha_i) = alpha_use;
	(*beta_i) = (*alpha_i) * multiplier;
	break;
      case 2:			/*should not happen */
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  scopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dge_sum_mv_testgen(int norm, enum blas_order_type order,
			     int m, int n, int randomize,
			     double *alpha, int alpha_flag, double *beta,
			     int beta_flag, double *a, int lda, double *b,
			     int ldb, double *x, int incx,
			     double *alpha_use_ptr, double *a_use,
			     double *b_use, int *seed, double *head_r_true,
			     double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dge_sum_mv{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) double*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) double*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) double*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem;
  double beta_zero_fake;
  double a_elem;
  double x_elem;
  double head_r_true_elem, tail_r_true_elem;
  double multiplier;
  double divider;
  double alpha_use;

  double *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = alpha_use_ptr;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  double *b_i = b;
  double *a_use_i = a_use;
  double *b_use_i = b_use;
  double *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if ((*alpha_i) == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if ((*beta_i) == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if ((*beta_i) == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem = xrand(seed);
      beta_i[0] = y_elem;
    }
    alpha_use = (*beta_i);
  } else {
    if (!alpha_flag) {
      y_elem = xrand(seed);
      alpha_i[0] = y_elem;
    }
    alpha_use = (*alpha_i);
  }
  /* put in return value */
  (*alpha_use_ptr_i) = alpha_use;

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*alpha_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*beta_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }




    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem = 0.0;
      BLAS_ddot_testgen(2 * n_i, 0, 2 * n_i, norm,
			blas_no_conj, &alpha_use, 1,
			&beta_zero_fake, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }


  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem = 0.0;
	BLAS_ddot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec, a_vec, seed,
			  &y_elem, &head_r_true_elem, &tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem;
	tail_r_true[ri] = tail_r_true_elem;
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  (*alpha_use_ptr_i) = alpha_use;


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {
    {
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	  case 3:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
    case 3:
      (*beta_i) = alpha_use;
      (*alpha_i) = (*beta_i) * multiplier;
      break;
    case 2:			/*should not happen */
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {
      {
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	    case 3:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
      case 3:
	(*alpha_i) = alpha_use;
	(*beta_i) = (*alpha_i) * multiplier;
	break;
      case 2:			/*should not happen */
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  dcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_cge_sum_mv_testgen(int norm, enum blas_order_type order,
			     int m, int n, int randomize,
			     void *alpha, int alpha_flag, void *beta,
			     int beta_flag, void *a, int lda, void *b,
			     int ldb, void *x, int incx, void *alpha_use_ptr,
			     void *a_use, void *b_use, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cge_sum_mv{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  float y_elem[2];
  float beta_zero_fake[2];
  float a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  float multiplier;
  float divider;
  float alpha_use[2];

  float *a_vec;
  float *x_vec;

  float *alpha_use_ptr_i = (float *) alpha_use_ptr;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;
  float *a_use_i = (float *) a_use;
  float *b_use_i = (float *) b_use;
  float *x_i = (float *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      float *a_vec_2;
      float *x_vec_2;
      a_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = xrand(seed);
	x_elem[1] = xrand(seed);
	x_elem[1] = 0.0;
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }
      /*copy new x into x_vec_2 (twice) */
      scopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      scopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_sdot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj,
			  alpha_use, 1,
			  beta_zero_fake, 1,
			  x_vec_2,
			  a_vec_2, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);


	/*multiply truth by 1+i (we will multiply 1+i to x later) */
	head_r_true_elem[1] = head_r_true_elem[0];
	tail_r_true_elem[1] = tail_r_true_elem[0];
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */

      /* also multiply x by 1+i, to compensate for change in
         truth above */
      for (j = 0; j < n_i; j++) {
	x_vec[2 * j] = x_vec_2[j];
	x_vec[2 * j + 1] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj, &alpha_use, 1,
			  &beta_zero_fake, 1, x_vec, a_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	b_i[aij] = a_elem[0];
	b_i[aij + 1] = a_elem[1];
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec, a_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      float one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      head_a_elem_2[0] = (double) a_elem[0] * divider;
	      tail_a_elem_2[0] = 0.0;
	      head_a_elem_2[1] = (double) a_elem[1] * divider;
	      tail_a_elem_2[1] = 0.0;
	    }
	    {
	      double cd[2];
	      cd[0] = (double) one_minus_i[0];
	      cd[1] = (double) one_minus_i[1];
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a0 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a1 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a1 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a0 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	    }
	    ((float *) a_elem)[0] = head_a_elem_3[0];
	    ((float *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	float one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		head_a_elem_2[0] = (double) a_elem[0] * divider;
		tail_a_elem_2[0] = 0.0;
		head_a_elem_2[1] = (double) a_elem[1] * divider;
		tail_a_elem_2[1] = 0.0;
	      }
	      {
		double cd[2];
		cd[0] = (double) one_minus_i[0];
		cd[1] = (double) one_minus_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_a_elem_2[0];
		  tail_a0 = tail_a_elem_2[0];
		  head_a1 = head_a_elem_2[1];
		  tail_a1 = tail_a_elem_2[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[0] = head_t1;
		  tail_a_elem_3[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[1] = head_t1;
		  tail_a_elem_3[1] = tail_t1;
		}

	      }
	      ((float *) a_elem)[0] = head_a_elem_3[0];
	      ((float *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  ccopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_testgen(int norm, enum blas_order_type order,
			     int m, int n, int randomize,
			     void *alpha, int alpha_flag, void *beta,
			     int beta_flag, void *a, int lda, void *b,
			     int ldb, void *x, int incx, void *alpha_use_ptr,
			     void *a_use, void *b_use, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  double a_elem[2];
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  double *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *b_i = (double *) b;
  double *a_use_i = (double *) a_use;
  double *b_use_i = (double *) b_use;
  double *x_i = (double *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_testgen(n_i, 0, n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      double *a_vec_2;
      double *x_vec_2;
      a_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = xrand(seed);
	x_elem[1] = xrand(seed);
	x_elem[1] = 0.0;
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }
      /*copy new x into x_vec_2 (twice) */
      dcopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      dcopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_ddot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj,
			  alpha_use, 1,
			  beta_zero_fake, 1,
			  x_vec_2,
			  a_vec_2, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);


	/*multiply truth by 1+i (we will multiply 1+i to x later) */
	head_r_true_elem[1] = head_r_true_elem[0];
	tail_r_true_elem[1] = tail_r_true_elem[0];
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */

      /* also multiply x by 1+i, to compensate for change in
         truth above */
      for (j = 0; j < n_i; j++) {
	x_vec[2 * j] = x_vec_2[j];
	x_vec[2 * j + 1] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj, &alpha_use, 1,
			  &beta_zero_fake, 1, x_vec, a_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    if (case_type == 3) {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    } else {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec, a_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      double one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = divider * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[0] = head_t;
	      tail_a_elem_2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = divider * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[1] = head_t;
	      tail_a_elem_2[1] = tail_t;
	    }
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_a_elem_2[0];
	      tail_a0 = tail_a_elem_2[0];
	      head_a1 = head_a_elem_2[1];
	      tail_a1 = tail_a_elem_2[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a0 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a1 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[0] = head_t1;
	      tail_a_elem_3[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a1 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a0 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[1] = head_t1;
	      tail_a_elem_3[1] = tail_t1;
	    }

	    ((double *) a_elem)[0] = head_a_elem_3[0];
	    ((double *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	double one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = divider * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[0] = head_t;
		tail_a_elem_2[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = divider * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[1] = head_t;
		tail_a_elem_2[1] = tail_t;
	      }
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a0 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a1 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a1 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a0 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	      ((double *) a_elem)[0] = head_a_elem_3[0];
	      ((double *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  zcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_cge_sum_mv_s_s_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, float *a, int lda, float *b,
				 int ldb, float *x, int incx,
				 void *alpha_use_ptr, float *a_use,
				 float *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cge_sum_mv_s_s{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) float*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) float*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    THIS CASE IS NOT PROPERLY TESTED.
 *      because of the difficulty in testing this case,
 *      a call with this case and randomize = 0 is
 *      converted into a call with randomize = 1.
 *      THERE IS INSUFFICIENT TESTING OF CANCELLATION IN THIS CASE.
 *      It is suggested that implementors be aware of this
 *      and take caution when working on ge_sum_mv.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  float y_elem[2];
  float beta_zero_fake[2];
  float a_elem;
  float x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];
  float multiplier;
  float divider;
  float alpha_use[2];

  float *a_vec;
  float *x_vec;

  float *alpha_use_ptr_i = (float *) alpha_use_ptr;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *b_i = b;
  float *a_use_i = a_use;
  float *b_use_i = b_use;
  float *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_s_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_s_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }

    /* case 3, Not fully tested */
    if (case_type == 3) {
      BLAS_cge_sum_mv_s_s_testgen(norm, order, m, n, 1 /*randomize */ ,
				  alpha_i, alpha_flag, beta_i, beta_flag,
				  a, lda, b, ldb, x, incx,
				  alpha_use_ptr_i, a_use, b_use,
				  seed, head_r_true, tail_r_true);
      blas_free(a_vec);
      blas_free(x_vec);
      return;
    }


    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem[0] = y_elem[1] = 0.0;
      BLAS_cdot_s_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


  } else {
    /* randomize == 1 */
    float *aa_vec;
    float *xx_vec;

    aa_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
    if (2 * n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    xx_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
    if (2 * n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < 2 * n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);
	blas_free(xx_vec);
	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);
	blas_free(xx_vec);
	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

	{
	  /* promote to complex */
	  int r;
	  for (r = 0; r < 2 * n_i; r++) {
	    aa_vec[2 * r] = a_vec[r];
	    aa_vec[2 * r + 1] = 0.0;
	  }
	}

	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  xx_vec,
			  aa_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }
    blas_free(aa_vec);
    blas_free(xx_vec);
  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 3:

	    a_elem = a_elem * divider;
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 3:

	      a_elem = a_elem * divider;
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  scopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_cge_sum_mv_s_c_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, float *a, int lda, float *b,
				 int ldb, void *x, int incx,
				 void *alpha_use_ptr, float *a_use,
				 float *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cge_sum_mv_s_c{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) float*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) float*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    THIS CASE IS NOT PROPERLY TESTED.
 *      because of the difficulty in testing this case,
 *      a call with this case and randomize = 0 is
 *      converted into a call with randomize = 1.
 *      THERE IS INSUFFICIENT TESTING OF CANCELLATION IN THIS CASE.
 *      It is suggested that implementors be aware of this
 *      and take caution when working on ge_sum_mv.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  float y_elem[2];
  float beta_zero_fake[2];
  float a_elem;
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  float multiplier;
  float divider;
  float alpha_use[2];

  float *a_vec;
  float *x_vec;

  float *alpha_use_ptr_i = (float *) alpha_use_ptr;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *b_i = b;
  float *a_use_i = a_use;
  float *b_use_i = b_use;
  float *x_i = (float *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_c_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_c_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }

    /* case 3, Not fully tested */
    if (case_type == 3) {
      BLAS_cge_sum_mv_s_c_testgen(norm, order, m, n, 1 /*randomize */ ,
				  alpha_i, alpha_flag, beta_i, beta_flag,
				  a, lda, b, ldb, x, incx,
				  alpha_use_ptr_i, a_use, b_use,
				  seed, head_r_true, tail_r_true);
      blas_free(a_vec);
      blas_free(x_vec);
      return;
    }


    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem[0] = y_elem[1] = 0.0;
      BLAS_cdot_c_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


  } else {
    /* randomize == 1 */
    float *aa_vec;


    aa_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
    if (2 * n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);

	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);

	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

	{
	  /* promote to complex */
	  int r;
	  for (r = 0; r < 2 * n_i; r++) {
	    aa_vec[2 * r] = a_vec[r];
	    aa_vec[2 * r + 1] = 0.0;
	  }
	}

	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec,
			  aa_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }
    blas_free(aa_vec);

  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 3:

	    a_elem = a_elem * divider;
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 3:

	      a_elem = a_elem * divider;
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  ccopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_cge_sum_mv_c_s_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, void *a, int lda, void *b,
				 int ldb, float *x, int incx,
				 void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cge_sum_mv_c_s{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  float y_elem[2];
  float beta_zero_fake[2];
  float a_elem[2];
  float x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];
  float multiplier;
  float divider;
  float alpha_use[2];

  float *a_vec;
  float *x_vec;

  float *alpha_use_ptr_i = (float *) alpha_use_ptr;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;
  float *a_use_i = (float *) a_use;
  float *b_use_i = (float *) b_use;
  float *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_s_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_cdot_s_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      float *a_vec_2;
      float *x_vec_2;
      a_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float));
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem = xrand(seed);

	x_i[xi] = x_elem;
      }
      /*copy new x into x_vec_2 (twice) */
      scopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      scopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_sdot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj,
			  alpha_use, 1,
			  beta_zero_fake, 1,
			  x_vec_2,
			  a_vec_2, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true_elem[1] = tail_r_true_elem[1] = 0.0;
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */
      for (j = 0; j < n_i; j++) {
	x_vec[j] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_s_c_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj, &alpha_use, 1,
			      &beta_zero_fake, 1, x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */

    float *xx_vec;


    xx_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
    if (2 * n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }
    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	b_i[aij] = a_elem[0];
	b_i[aij + 1] = a_elem[1];
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < 2 * n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    a_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);

	blas_free(xx_vec);
	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    a_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);

	blas_free(xx_vec);
	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_cdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  xx_vec,
			  a_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(xx_vec);
  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      float one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      head_a_elem_2[0] = (double) a_elem[0] * divider;
	      tail_a_elem_2[0] = 0.0;
	      head_a_elem_2[1] = (double) a_elem[1] * divider;
	      tail_a_elem_2[1] = 0.0;
	    }
	    {
	      double cd[2];
	      cd[0] = (double) one_minus_i[0];
	      cd[1] = (double) one_minus_i[1];
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a0 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a1 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a1 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a0 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	    }
	    ((float *) a_elem)[0] = head_a_elem_3[0];
	    ((float *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	float one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		head_a_elem_2[0] = (double) a_elem[0] * divider;
		tail_a_elem_2[0] = 0.0;
		head_a_elem_2[1] = (double) a_elem[1] * divider;
		tail_a_elem_2[1] = 0.0;
	      }
	      {
		double cd[2];
		cd[0] = (double) one_minus_i[0];
		cd[1] = (double) one_minus_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_a_elem_2[0];
		  tail_a0 = tail_a_elem_2[0];
		  head_a1 = head_a_elem_2[1];
		  tail_a1 = tail_a_elem_2[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[0] = head_t1;
		  tail_a_elem_3[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[1] = head_t1;
		  tail_a_elem_3[1] = tail_t1;
		}

	      }
	      ((float *) a_elem)[0] = head_a_elem_3[0];
	      ((float *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  scopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_d_d_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, double *a, int lda, double *b,
				 int ldb, double *x, int incx,
				 void *alpha_use_ptr, double *a_use,
				 double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_d_d{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) double*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) double*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    THIS CASE IS NOT PROPERLY TESTED.
 *      because of the difficulty in testing this case,
 *      a call with this case and randomize = 0 is
 *      converted into a call with randomize = 1.
 *      THERE IS INSUFFICIENT TESTING OF CANCELLATION IN THIS CASE.
 *      It is suggested that implementors be aware of this
 *      and take caution when working on ge_sum_mv.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  double a_elem;
  double x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  double *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *b_i = b;
  double *a_use_i = a_use;
  double *b_use_i = b_use;
  double *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_d_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_d_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }

    /* case 3, Not fully tested */
    if (case_type == 3) {
      BLAS_zge_sum_mv_d_d_testgen(norm, order, m, n, 1 /*randomize */ ,
				  alpha_i, alpha_flag, beta_i, beta_flag,
				  a, lda, b, ldb, x, incx,
				  alpha_use_ptr_i, a_use, b_use,
				  seed, head_r_true, tail_r_true);
      blas_free(a_vec);
      blas_free(x_vec);
      return;
    }


    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem[0] = y_elem[1] = 0.0;
      BLAS_zdot_d_d_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


  } else {
    /* randomize == 1 */
    double *aa_vec;
    double *xx_vec;

    aa_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
    if (2 * n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    xx_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
    if (2 * n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < 2 * n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);
	blas_free(xx_vec);
	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);
	blas_free(xx_vec);
	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

	{
	  /* promote to complex */
	  int r;
	  for (r = 0; r < 2 * n_i; r++) {
	    aa_vec[2 * r] = a_vec[r];
	    aa_vec[2 * r + 1] = 0.0;
	  }
	}

	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  xx_vec,
			  aa_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }
    blas_free(aa_vec);
    blas_free(xx_vec);
  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 3:

	    a_elem = a_elem * divider;
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 3:

	      a_elem = a_elem * divider;
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  dcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_d_z_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, double *a, int lda, double *b,
				 int ldb, void *x, int incx,
				 void *alpha_use_ptr, double *a_use,
				 double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_d_z{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) double*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) double*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    THIS CASE IS NOT PROPERLY TESTED.
 *      because of the difficulty in testing this case,
 *      a call with this case and randomize = 0 is
 *      converted into a call with randomize = 1.
 *      THERE IS INSUFFICIENT TESTING OF CANCELLATION IN THIS CASE.
 *      It is suggested that implementors be aware of this
 *      and take caution when working on ge_sum_mv.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  double a_elem;
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  double *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *b_i = b;
  double *a_use_i = a_use;
  double *b_use_i = b_use;
  double *x_i = (double *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_z_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_z_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }

    /* case 3, Not fully tested */
    if (case_type == 3) {
      BLAS_zge_sum_mv_d_z_testgen(norm, order, m, n, 1 /*randomize */ ,
				  alpha_i, alpha_flag, beta_i, beta_flag,
				  a, lda, b, ldb, x, incx,
				  alpha_use_ptr_i, a_use, b_use,
				  seed, head_r_true, tail_r_true);
      blas_free(a_vec);
      blas_free(x_vec);
      return;
    }


    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem[0] = y_elem[1] = 0.0;
      BLAS_zdot_z_d_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


  } else {
    /* randomize == 1 */
    double *aa_vec;


    aa_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
    if (2 * n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);

	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  {
	    /* promote to complex */
	    int r;
	    for (r = 0; r < n_i; r++) {
	      aa_vec[2 * r] = a_vec[r];
	      aa_vec[2 * r + 1] = 0.0;
	    }
	  }
	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    x_vec,
			    aa_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	blas_free(aa_vec);

	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

	{
	  /* promote to complex */
	  int r;
	  for (r = 0; r < 2 * n_i; r++) {
	    aa_vec[2 * r] = a_vec[r];
	    aa_vec[2 * r + 1] = 0.0;
	  }
	}

	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  x_vec,
			  aa_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }
    blas_free(aa_vec);

  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 3:

	    a_elem = a_elem * divider;
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 3:

	      a_elem = a_elem * divider;
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  zcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_z_d_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, void *a, int lda, void *b,
				 int ldb, double *x, int incx,
				 void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_z_d{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  double a_elem[2];
  double x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  double *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *b_i = (double *) b;
  double *a_use_i = (double *) a_use;
  double *b_use_i = (double *) b_use;
  double *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_d_z_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_d_z_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      double *a_vec_2;
      double *x_vec_2;
      a_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double));
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem = xrand(seed);

	x_i[xi] = x_elem;
      }
      /*copy new x into x_vec_2 (twice) */
      dcopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      dcopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_ddot_testgen(2 * n_i, 0, 2 * n_i, norm,
			  blas_no_conj,
			  alpha_use, 1,
			  beta_zero_fake, 1,
			  x_vec_2,
			  a_vec_2, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true_elem[1] = tail_r_true_elem[1] = 0.0;
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */
      for (j = 0; j < n_i; j++) {
	x_vec[j] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_d_z_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj, &alpha_use, 1,
			      &beta_zero_fake, 1, x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */

    double *xx_vec;


    xx_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
    if (2 * n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    if (case_type == 3) {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    } else {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < 2 * n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    a_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);

	blas_free(xx_vec);
	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj,
			    &alpha_use, 1,
			    &beta_zero_fake, 1,
			    xx_vec,
			    a_vec,
			    seed, y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);

	blas_free(xx_vec);
	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			  &alpha_use, 1,
			  &beta_zero_fake, 1,
			  xx_vec,
			  a_vec,
			  seed, y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }

    blas_free(xx_vec);
  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      double one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = divider * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[0] = head_t;
	      tail_a_elem_2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = divider * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[1] = head_t;
	      tail_a_elem_2[1] = tail_t;
	    }
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_a_elem_2[0];
	      tail_a0 = tail_a_elem_2[0];
	      head_a1 = head_a_elem_2[1];
	      tail_a1 = tail_a_elem_2[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a0 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a1 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[0] = head_t1;
	      tail_a_elem_3[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a1 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a0 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[1] = head_t1;
	      tail_a_elem_3[1] = tail_t1;
	    }

	    ((double *) a_elem)[0] = head_a_elem_3[0];
	    ((double *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	double one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = divider * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[0] = head_t;
		tail_a_elem_2[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = divider * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[1] = head_t;
		tail_a_elem_2[1] = tail_t;
	      }
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a0 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a1 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a1 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a0 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	      ((double *) a_elem)[0] = head_a_elem_3[0];
	      ((double *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  dcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dge_sum_mv_s_s_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 double *alpha, int alpha_flag, double *beta,
				 int beta_flag, float *a, int lda, float *b,
				 int ldb, float *x, int incx,
				 double *alpha_use_ptr, float *a_use,
				 float *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dge_sum_mv_s_s{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) double*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) float*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) float*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem;
  double beta_zero_fake;
  float a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;
  double multiplier;
  double divider;
  double alpha_use;

  float *a_vec;
  float *x_vec;

  double *alpha_use_ptr_i = alpha_use_ptr;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  float *b_i = b;
  float *a_use_i = a_use;
  float *b_use_i = b_use;
  float *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if ((*alpha_i) == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if ((*beta_i) == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if ((*beta_i) == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }
    alpha_use = (*beta_i);
  } else {
    if (!alpha_flag) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    alpha_use = (*alpha_i);
  }
  /* put in return value */
  (*alpha_use_ptr_i) = alpha_use;

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_s_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*alpha_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_s_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*beta_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }




    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem = 0.0;
      BLAS_ddot_s_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }


  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_s_s_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_s_s_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem = 0.0;
	BLAS_ddot_s_s_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      &y_elem, &head_r_true_elem, &tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem;
	tail_r_true[ri] = tail_r_true_elem;
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  (*alpha_use_ptr_i) = alpha_use;


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {
    {
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	  case 3:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
    case 3:
      (*beta_i) = alpha_use;
      (*alpha_i) = (*beta_i) * multiplier;
      break;
    case 2:			/*should not happen */
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {
      {
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	    case 3:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
      case 3:
	(*alpha_i) = alpha_use;
	(*beta_i) = (*alpha_i) * multiplier;
	break;
      case 2:			/*should not happen */
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  scopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dge_sum_mv_s_d_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 double *alpha, int alpha_flag, double *beta,
				 int beta_flag, float *a, int lda, float *b,
				 int ldb, double *x, int incx,
				 double *alpha_use_ptr, float *a_use,
				 float *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dge_sum_mv_s_d{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) double*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) float*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) float*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem;
  double beta_zero_fake;
  float a_elem;
  double x_elem;
  double head_r_true_elem, tail_r_true_elem;
  double multiplier;
  double divider;
  double alpha_use;

  float *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = alpha_use_ptr;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  float *b_i = b;
  float *a_use_i = a_use;
  float *b_use_i = b_use;
  double *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if ((*alpha_i) == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if ((*beta_i) == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if ((*beta_i) == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }
    alpha_use = (*beta_i);
  } else {
    if (!alpha_flag) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    alpha_use = (*alpha_i);
  }
  /* put in return value */
  (*alpha_use_ptr_i) = alpha_use;

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_d_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*alpha_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_d_s_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*beta_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }




    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem = 0.0;
      BLAS_ddot_d_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      sge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      sge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }


  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    dcopy_vector(x, n_i, incx, x_vec, 1);
    dcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_d_s_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_d_s_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	dcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	sge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	sge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem = 0.0;
	BLAS_ddot_d_s_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      &y_elem, &head_r_true_elem, &tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem;
	tail_r_true[ri] = tail_r_true_elem;
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  (*alpha_use_ptr_i) = alpha_use;


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {
    {
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	  case 3:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
    case 3:
      (*beta_i) = alpha_use;
      (*alpha_i) = (*beta_i) * multiplier;
      break;
    case 2:			/*should not happen */
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {
      {
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	    case 3:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
      case 3:
	(*alpha_i) = alpha_use;
	(*beta_i) = (*alpha_i) * multiplier;
	break;
      case 2:			/*should not happen */
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  dcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dge_sum_mv_d_s_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 double *alpha, int alpha_flag, double *beta,
				 int beta_flag, double *a, int lda, double *b,
				 int ldb, float *x, int incx,
				 double *alpha_use_ptr, double *a_use,
				 double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dge_sum_mv_d_s{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) double*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) double*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) double*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem;
  double beta_zero_fake;
  double a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;
  double multiplier;
  double divider;
  double alpha_use;

  double *a_vec;
  float *x_vec;

  double *alpha_use_ptr_i = alpha_use_ptr;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  double *b_i = b;
  double *a_use_i = a_use;
  double *b_use_i = b_use;
  float *x_i = x;

  n_i = n;
  m_i = m;

  beta_zero_fake = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;


  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if ((*alpha_i) == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if ((*beta_i) == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if ((*beta_i) == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }
    alpha_use = (*beta_i);
  } else {
    if (!alpha_flag) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    alpha_use = (*alpha_i);
  }
  /* put in return value */
  (*alpha_use_ptr_i) = alpha_use;

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_s_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    b_i[aij] = a_elem;
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*alpha_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem = 0.0;
	  BLAS_ddot_s_d_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = (float) xrand(seed);
	    a_i[aij] = a_elem;
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	(*beta_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }




    /* Fill in matrix A, B */
    for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
      y_elem = 0.0;
      BLAS_ddot_s_d_testgen(2 * n_i, 0, 2 * n_i, norm,
			    blas_no_conj, &alpha_use, 1,
			    &beta_zero_fake, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      dge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      dge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);

      /*commits an element to the truth */
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }


  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem = (float) xrand(seed);
      x_i[xi * incxi] = x_elem;
    }


    /*set a randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    /*set b randomly */
    if (order == blas_colmajor) {
      incai = 1;
      incaij = ldb;
    } else {
      incai = ldb;
      incaij = 1;
    }



    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	b_i[aij] = a_elem;
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    scopy_vector(x, n_i, incx, x_vec, 1);
    scopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_s_d_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    b_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem = 0.0;

	  BLAS_ddot_s_d_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				&y_elem,
				&head_r_true_elem, &tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem;
	  tail_r_true[ri] = tail_r_true_elem;
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = a_i[aij];
	    a_use_i[aij] = a_elem;
	  }
	}
	(*alpha_use_ptr_i) = alpha_use;
	scopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	dge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	dge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem = 0.0;
	BLAS_ddot_s_d_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      &y_elem, &head_r_true_elem, &tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem;
	tail_r_true[ri] = tail_r_true_elem;
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = a_i[aij];
      a_use_i[aij] = a_elem;
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }



  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem = b_i[aij];
      b_use_i[aij] = a_elem;
    }
  }
  (*alpha_use_ptr_i) = alpha_use;


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {
    {
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }



      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem = a_i[aij];
	  switch (case_type) {
	  case 1:
	  case 3:
	    a_elem = a_elem * divider;
	    break;
	  case 2:		/*should not happen */
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem;
	}
      }
    }

    switch (case_type) {
    case 1:
    case 3:
      (*beta_i) = alpha_use;
      (*alpha_i) = (*beta_i) * multiplier;
      break;
    case 2:			/*should not happen */
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {
      {
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}



	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem = b_i[aij];
	    switch (case_type) {
	    case 1:
	    case 3:
	      a_elem = a_elem * divider;
	      break;
	    case 2:		/*should not happen */
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem;
	  }
	}
      }

      switch (case_type) {
      case 1:
      case 3:
	(*alpha_i) = alpha_use;
	(*beta_i) = (*alpha_i) * multiplier;
	break;
      case 2:			/*should not happen */
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  scopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_c_c_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, void *a, int lda, void *b,
				 int ldb, void *x, int incx,
				 void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_c_c{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  float a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  float *a_vec;
  float *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;
  float *a_use_i = (float *) a_use;
  float *b_use_i = (float *) b_use;
  float *x_i = (float *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_c_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_c_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      float *a_vec_2;
      float *x_vec_2;
      a_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_elem[1] = 0.0;
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }
      /*copy new x into x_vec_2 (twice) */
      scopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      scopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_ddot_s_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj,
			      alpha_use, 1,
			      beta_zero_fake, 1,
			      x_vec_2,
			      a_vec_2, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);


	/*multiply truth by 1+i (we will multiply 1+i to x later) */
	head_r_true_elem[1] = head_r_true_elem[0];
	tail_r_true_elem[1] = tail_r_true_elem[0];
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */

      /* also multiply x by 1+i, to compensate for change in
         truth above */
      for (j = 0; j < n_i; j++) {
	x_vec[2 * j] = x_vec_2[j];
	x_vec[2 * j + 1] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_c_c_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj, &alpha_use, 1,
			      &beta_zero_fake, 1, x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    if (case_type == 3) {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    } else {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_c_c_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_c_c_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_c_c_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      float one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      head_a_elem_2[0] = (double) a_elem[0] * divider;
	      tail_a_elem_2[0] = 0.0;
	      head_a_elem_2[1] = (double) a_elem[1] * divider;
	      tail_a_elem_2[1] = 0.0;
	    }
	    {
	      double cd[2];
	      cd[0] = (double) one_minus_i[0];
	      cd[1] = (double) one_minus_i[1];
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a0 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a1 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a1 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a0 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	    }
	    ((float *) a_elem)[0] = head_a_elem_3[0];
	    ((float *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	float one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		head_a_elem_2[0] = (double) a_elem[0] * divider;
		tail_a_elem_2[0] = 0.0;
		head_a_elem_2[1] = (double) a_elem[1] * divider;
		tail_a_elem_2[1] = 0.0;
	      }
	      {
		double cd[2];
		cd[0] = (double) one_minus_i[0];
		cd[1] = (double) one_minus_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_a_elem_2[0];
		  tail_a0 = tail_a_elem_2[0];
		  head_a1 = head_a_elem_2[1];
		  tail_a1 = tail_a_elem_2[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[0] = head_t1;
		  tail_a_elem_3[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[1] = head_t1;
		  tail_a_elem_3[1] = tail_t1;
		}

	      }
	      ((float *) a_elem)[0] = head_a_elem_3[0];
	      ((float *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  ccopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_c_z_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, void *a, int lda, void *b,
				 int ldb, void *x, int incx,
				 void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_c_z{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  float a_elem[2];
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  float *a_vec;
  double *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  float *b_i = (float *) b;
  float *a_use_i = (float *) a_use;
  float *b_use_i = (float *) b_use;
  double *x_i = (double *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_z_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_z_c_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      float *a_vec_2;
      double *x_vec_2;
      a_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_elem[1] = 0.0;
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }
      /*copy new x into x_vec_2 (twice) */
      dcopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      dcopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_ddot_d_s_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj,
			      alpha_use, 1,
			      beta_zero_fake, 1,
			      x_vec_2,
			      a_vec_2, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);


	/*multiply truth by 1+i (we will multiply 1+i to x later) */
	head_r_true_elem[1] = head_r_true_elem[0];
	tail_r_true_elem[1] = tail_r_true_elem[0];
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */

      /* also multiply x by 1+i, to compensate for change in
         truth above */
      for (j = 0; j < n_i; j++) {
	x_vec[2 * j] = x_vec_2[j];
	x_vec[2 * j + 1] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_z_c_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj, &alpha_use, 1,
			      &beta_zero_fake, 1, x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	cge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    if (case_type == 3) {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    } else {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    zcopy_vector(x, n_i, incx, x_vec, 1);
    zcopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_z_c_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_z_c_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	zcopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	cge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	cge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_z_c_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      float one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      head_a_elem_2[0] = (double) a_elem[0] * divider;
	      tail_a_elem_2[0] = 0.0;
	      head_a_elem_2[1] = (double) a_elem[1] * divider;
	      tail_a_elem_2[1] = 0.0;
	    }
	    {
	      double cd[2];
	      cd[0] = (double) one_minus_i[0];
	      cd[1] = (double) one_minus_i[1];
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a0 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a1 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = cd[0] * split;
		  b1 = con - cd[0];
		  b1 = con - b1;
		  b2 = cd[0] - b1;

		  c11 = head_a1 * cd[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * cd[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = cd[1] * split;
		  b1 = con - cd[1];
		  b1 = con - b1;
		  b2 = cd[1] - b1;

		  c11 = head_a0 * cd[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * cd[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	    }
	    ((float *) a_elem)[0] = head_a_elem_3[0];
	    ((float *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	float one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		head_a_elem_2[0] = (double) a_elem[0] * divider;
		tail_a_elem_2[0] = 0.0;
		head_a_elem_2[1] = (double) a_elem[1] * divider;
		tail_a_elem_2[1] = 0.0;
	      }
	      {
		double cd[2];
		cd[0] = (double) one_minus_i[0];
		cd[1] = (double) one_minus_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_a_elem_2[0];
		  tail_a0 = tail_a_elem_2[0];
		  head_a1 = head_a_elem_2[1];
		  tail_a1 = tail_a_elem_2[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[0] = head_t1;
		  tail_a_elem_3[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_a_elem_3[1] = head_t1;
		  tail_a_elem_3[1] = tail_t1;
		}

	      }
	      ((float *) a_elem)[0] = head_a_elem_3[0];
	      ((float *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  zcopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zge_sum_mv_z_c_testgen(int norm, enum blas_order_type order,
				 int m, int n, int randomize,
				 void *alpha, int alpha_flag, void *beta,
				 int beta_flag, void *a, int lda, void *b,
				 int ldb, void *x, int incx,
				 void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zge_sum_mv_z_c{_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
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
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) void*
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) void*
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) void*
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
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
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * 
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 */
{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  double y_elem[2];
  double beta_zero_fake[2];
  double a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double multiplier;
  double divider;
  double alpha_use[2];

  double *a_vec;
  float *x_vec;

  double *alpha_use_ptr_i = (double *) alpha_use_ptr;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *b_i = (double *) b;
  double *a_use_i = (double *) a_use;
  double *b_use_i = (double *) b_use;
  float *x_i = (float *) x;

  n_i = n;
  m_i = m;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  /*x_vec, a_vec must have stride of 1 */
  inca_veci = 1;
  inca_veci *= 2;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  incri = 1;
  incri *= 2;

  incxi = incx;
  incx_veci = 1;
  incx_veci *= 2;
  incxi *= 2;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  /* choose k */
  k = 0;
  while (!k) {
    k = xrand(seed) * 7 - 3;
  }

  multiplier = 1.0;
  divider = 1.0;
  for (i = 0; i < k; i++) {
    multiplier = multiplier * 2.0;
    divider = divider * 0.5;
  }
  for (i = 0; i > k; i--) {
    multiplier = multiplier * 0.5;
    divider = divider * 2.0;
  }
  /* decide which case */
  if (alpha_flag) {
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* case 2 */
      case_type = 2;
      which_free = ALPHA_USE_IS_BETA;	/* for use beta */
    } else {
      if (beta_flag) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* case 2 */
	  case_type = 2;
	  which_free = ALPHA_USE_IS_ALPHA;
	  /*for use alpha */
	} else {
	  /* case 4 */
	  case_type = 4;
	  k = 0;
	  which_free = ALPHA_USE_IS_EITHER;
	}
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_ALPHA;
	/* for beta free, use alpha */
      }
    }
  } else {
    if (beta_flag) {
      if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	/* case 2 */
	case_type = 2;
	which_free = ALPHA_USE_IS_ALPHA;
	/*alpha is nonzero */
      } else {
	/* case 3 */
	case_type = 3;
	which_free = ALPHA_USE_IS_BETA;
	/* for alpha free, use beta */
      }
    } else {
      /* case 1 */
      case_type = 1;
      which_free = ALPHA_USE_IS_ALPHA;
    }
  }

  if (which_free == ALPHA_USE_IS_BETA) {
    if (!beta_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = beta_i[0];
    alpha_use[1] = beta_i[1];
  } else {
    if (!alpha_flag) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    alpha_use[0] = alpha_i[0];
    alpha_use[1] = alpha_i[1];
  }
  /* put in return value */
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];

  if (randomize == 0) {

    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);

    if (case_type == 2) {
      /* degenerate case - similar to gemv */
      if (which_free == ALPHA_USE_IS_ALPHA) {
	/* alpha == alpha_use */

	/* now Fill in matrix alpha only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_c_z_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill a, x, and return */

	/*set b randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      } else {

	/* now Fill in matrix beta only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {

	  y_elem[0] = y_elem[1] = 0.0;
	  BLAS_zdot_c_z_testgen(n_i, 0, n_i, norm,
				blas_no_conj, &alpha_use, 1,
				&beta_zero_fake, 1, x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);

	  /*commits an element to the truth */
	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*now fill b, x, and return */

	/*set a randomly */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = (float) xrand(seed);
	    a_elem[1] = (float) xrand(seed);
	    a_i[aij] = a_elem[0];
	    a_i[aij + 1] = a_elem[1];
	  }
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	beta_i[0] = alpha_use[0];
	beta_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);
	return;
      }
    }



    /* case 3, start with real matricies, x */
    if (case_type == 3) {
      double *a_vec_2;
      float *x_vec_2;
      a_vec_2 = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
      if (4 * n_i > 0 && a_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      x_vec_2 = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
      if (4 * n_i > 0 && x_vec_2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < 2 * n_i * inca_veci; i += inca_veci) {
	a_vec[i] = 0.0;
	a_vec[i + 1] = 0.0;
      }

      /*first pick x randomly, but real */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_elem[1] = 0.0;
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }
      /*copy new x into x_vec_2 (twice) */
      scopy_vector(x, n_i, 2 * incx, x_vec_2, 1);
      scopy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

      /* Now Fill in matrix A, B real */
      /*since we have case 3, we know alpha_use == 1.0+0i,
         so we will force it to be real */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_ddot_s_d_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj,
			      alpha_use, 1,
			      beta_zero_fake, 1,
			      x_vec_2,
			      a_vec_2, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);


	/*multiply truth by 1+i (we will multiply 1+i to x later) */
	head_r_true_elem[1] = head_r_true_elem[0];
	tail_r_true_elem[1] = tail_r_true_elem[0];
	for (j = 0; j < 2 * n_i; j++) {
	  a_vec[2 * j] = a_vec_2[j];
	  a_vec[2 * j + 1] = 0.0;
	}
	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       a_vec + inca_veci * n_i, i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
      /* copy to x_vec - will be copied to x_i later */

      /* also multiply x by 1+i, to compensate for change in
         truth above */
      for (j = 0; j < n_i; j++) {
	x_vec[2 * j] = x_vec_2[j];
	x_vec[2 * j + 1] = x_vec_2[j];
      }
      blas_free(x_vec_2);
      blas_free(a_vec_2);
    } else {
      /*not case 3 */

      /* Fill in matrix A, B */
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_c_z_testgen(2 * n_i, 0, 2 * n_i, norm,
			      blas_no_conj, &alpha_use, 1,
			      &beta_zero_fake, 1, x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	zge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb,
		       (a_vec + inca_veci * n_i), i);

	/*commits an element to the truth */
	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }

    }

  } else {
    /* randomize == 1 */







    /*first pick x randomly */
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi * incxi] = x_elem[0];
      x_i[xi * incxi + 1] = x_elem[1];
    }
    if (case_type == 3) {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = 0.0;
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    } else {

      /*set a randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }

      /*set b randomly */
      if (order == blas_colmajor) {
	incai = 1;
	incaij = ldb;
      } else {
	incai = ldb;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  b_i[aij] = a_elem[0];
	  b_i[aij + 1] = a_elem[1];
	}
      }
    }

    /* now compute appropriate truth */

    /* get x */
    /*copy new x into x_vec (twice) */
    ccopy_vector(x, n_i, incx, x_vec, 1);
    ccopy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);


    if (case_type == 2) {
      if (which_free == ALPHA_USE_IS_BETA) {

	/* Fill in truth from b, beta_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_c_z_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set a_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_use_i[aij] = 0.0;
	    a_use_i[aij + 1] = 0.0;
	  }
	}

	/*set b_use = b */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    b_use_i[aij] = a_elem[0];
	    b_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;

      } else {

	/* Fill in truth from a, alpha_i only */
	for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	  zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);


	  y_elem[0] = y_elem[1] = 0.0;

	  BLAS_zdot_c_z_testgen(n_i, n_i, 0, norm, blas_no_conj,
				&alpha_use, 1,
				&beta_zero_fake, 1,
				x_vec, a_vec, seed,
				y_elem, head_r_true_elem, tail_r_true_elem);

	  head_r_true[ri] = head_r_true_elem[0];
	  head_r_true[ri + 1] = head_r_true_elem[1];
	  tail_r_true[ri] = tail_r_true_elem[0];
	  tail_r_true[ri + 1] = tail_r_true_elem[1];
	}

	/*set b_use = 0 */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    b_use_i[aij] = 0.0;
	    b_use_i[aij + 1] = 0.0;
	  }
	}

	/*set a_use = a */
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = lda;
	} else {
	  incai = lda;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = a_i[aij];
	    a_elem[1] = a_i[aij + 1];
	    a_use_i[aij] = a_elem[0];
	    a_use_i[aij + 1] = a_elem[1];
	  }
	}
	alpha_use_ptr_i[0] = alpha_use[0];
	alpha_use_ptr_i[1] = alpha_use[1];
	ccopy_vector(x_vec, n_i, 1, x, incx);
	blas_free(a_vec);
	blas_free(x_vec);


	return;
      }
    } else {
      for (i = 0, ri = 0; i < m_i; i++, ri += incri) {
	zge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
	zge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb,
		     (a_vec + inca_veci * n_i), i);



	y_elem[0] = y_elem[1] = 0.0;
	BLAS_zdot_c_z_testgen(2 * n_i, 2 * n_i, 0, norm, blas_no_conj,
			      &alpha_use, 1,
			      &beta_zero_fake, 1,
			      x_vec, a_vec, seed,
			      y_elem, head_r_true_elem, tail_r_true_elem);

	head_r_true[ri] = head_r_true_elem[0];
	head_r_true[ri + 1] = head_r_true_elem[1];
	tail_r_true[ri] = tail_r_true_elem[0];
	tail_r_true[ri + 1] = tail_r_true_elem[1];
      }
    }


  }



  /*set a_use = a */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      a_use_i[aij] = a_elem[0];
      a_use_i[aij + 1] = a_elem[1];
    }
  }

  /*set b_use = b */
  if (order == blas_colmajor) {
    incai = 1;
    incaij = ldb;
  } else {
    incai = ldb;
    incaij = 1;
  }
  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
    for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
      a_elem[0] = b_i[aij];
      a_elem[1] = b_i[aij + 1];
      b_use_i[aij] = a_elem[0];
      b_use_i[aij + 1] = a_elem[1];
    }
  }
  alpha_use_ptr_i[0] = alpha_use[0];
  alpha_use_ptr_i[1] = alpha_use[1];


  /* now we scale */
  if (which_free == ALPHA_USE_IS_BETA) {

    {

      double one_minus_i[2];
      double head_a_elem_2[2], tail_a_elem_2[2];
      double head_a_elem_3[2], tail_a_elem_3[2];
      one_minus_i[0] = 0.5;
      one_minus_i[1] = -0.5;
      if (order == blas_colmajor) {
	incai = 1;
	incaij = lda;
      } else {
	incai = lda;
	incaij = 1;
      }
      incai *= 2;
      incaij *= 2;

      for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = a_i[aij];
	  a_elem[1] = a_i[aij + 1];
	  switch (case_type) {
	  case 1:
	    {
	      a_elem[0] = a_elem[0] * divider;
	      a_elem[1] = a_elem[1] * divider;
	    }
	    break;
	  case 2:		/*should not happen */
	  case 3:
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[0] * split;
		b1 = con - a_elem[0];
		b1 = con - b1;
		b2 = a_elem[0] - b1;

		head_t = divider * a_elem[0];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[0] = head_t;
	      tail_a_elem_2[0] = tail_t;
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = divider * split;
		a1 = con - divider;
		a1 = con - a1;
		a2 = divider - a1;
		con = a_elem[1] * split;
		b1 = con - a_elem[1];
		b1 = con - b1;
		b2 = a_elem[1] - b1;

		head_t = divider * a_elem[1];
		tail_t = (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_a_elem_2[1] = head_t;
	      tail_a_elem_2[1] = tail_t;
	    }
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_a_elem_2[0];
	      tail_a0 = tail_a_elem_2[0];
	      head_a1 = head_a_elem_2[1];
	      tail_a1 = tail_a_elem_2[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a0 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a1 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[0] = head_t1;
	      tail_a_elem_3[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = one_minus_i[0] * split;
		b1 = con - one_minus_i[0];
		b1 = con - b1;
		b2 = one_minus_i[0] - b1;

		c11 = head_a1 * one_minus_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * one_minus_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = one_minus_i[1] * split;
		b1 = con - one_minus_i[1];
		b1 = con - b1;
		b2 = one_minus_i[1] - b1;

		c11 = head_a0 * one_minus_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * one_minus_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_a_elem_3[1] = head_t1;
	      tail_a_elem_3[1] = tail_t1;
	    }

	    ((double *) a_elem)[0] = head_a_elem_3[0];
	    ((double *) a_elem)[1] = head_a_elem_3[1];
	    break;
	  case 4:		/*k ==0 */
	    break;
	  }
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	}
      }
    }

    switch (case_type) {
    case 1:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      break;
    case 2:			/*should not happen */
      break;
    case 3:
      beta_i[0] = alpha_use[0];
      beta_i[1] = alpha_use[1];
      {
	alpha_i[0] = beta_i[0] * multiplier;
	alpha_i[1] = beta_i[1] * multiplier;
      }
      alpha_i[1] = alpha_i[0];
      break;
    case 4:
      break;
    }
  } else {
    if (which_free == ALPHA_USE_IS_ALPHA) {

      {

	double one_minus_i[2];
	double head_a_elem_2[2], tail_a_elem_2[2];
	double head_a_elem_3[2], tail_a_elem_3[2];
	one_minus_i[0] = 0.5;
	one_minus_i[1] = -0.5;
	if (order == blas_colmajor) {
	  incai = 1;
	  incaij = ldb;
	} else {
	  incai = ldb;
	  incaij = 1;
	}
	incai *= 2;
	incaij *= 2;

	for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
	  for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	    a_elem[0] = b_i[aij];
	    a_elem[1] = b_i[aij + 1];
	    switch (case_type) {
	    case 1:
	      {
		a_elem[0] = a_elem[0] * divider;
		a_elem[1] = a_elem[1] * divider;
	      }
	      break;
	    case 2:		/*should not happen */
	    case 3:
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = divider * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[0] = head_t;
		tail_a_elem_2[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = divider * split;
		  a1 = con - divider;
		  a1 = con - a1;
		  a2 = divider - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = divider * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_a_elem_2[1] = head_t;
		tail_a_elem_2[1] = tail_t;
	      }
	      {
		/* Compute complex-extra = complex-extra * complex-double. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		head_a0 = head_a_elem_2[0];
		tail_a0 = tail_a_elem_2[0];
		head_a1 = head_a_elem_2[1];
		tail_a1 = tail_a_elem_2[1];
		/* real part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a0 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a1 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[0] = head_t1;
		tail_a_elem_3[0] = tail_t1;
		/* imaginary part */
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = one_minus_i[0] * split;
		  b1 = con - one_minus_i[0];
		  b1 = con - b1;
		  b2 = one_minus_i[0] - b1;

		  c11 = head_a1 * one_minus_i[0];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * one_minus_i[0];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = one_minus_i[1] * split;
		  b1 = con - one_minus_i[1];
		  b1 = con - b1;
		  b2 = one_minus_i[1] - b1;

		  c11 = head_a0 * one_minus_i[1];
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * one_minus_i[1];
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t2 = t1 + t2;
		  tail_t2 = t2 - (head_t2 - t1);
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_a_elem_3[1] = head_t1;
		tail_a_elem_3[1] = tail_t1;
	      }

	      ((double *) a_elem)[0] = head_a_elem_3[0];
	      ((double *) a_elem)[1] = head_a_elem_3[1];
	      break;
	    case 4:		/*k ==0 */
	      break;
	    }
	    b_i[aij] = a_elem[0];
	    b_i[aij + 1] = a_elem[1];
	  }
	}
      }

      switch (case_type) {
      case 1:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	break;
      case 2:			/*should not happen */
	break;
      case 3:
	alpha_i[0] = alpha_use[0];
	alpha_i[1] = alpha_use[1];
	{
	  beta_i[0] = alpha_i[0] * multiplier;
	  beta_i[1] = alpha_i[1] * multiplier;
	}
	beta_i[1] = beta_i[0];
	break;
      case 4:
	break;
      }
    } else {
      /*which_free = ALPHA_USE_IS_EITHER , case 4 */
    }
  }				/* which_free if */

  /*copy x_vec into x : it is possible that the generator
     changed x_vec, even though none were free */
  ccopy_vector(x_vec, n_i, 1, x, incx);
  blas_free(a_vec);
  blas_free(x_vec);
}
