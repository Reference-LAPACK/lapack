#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void BLAS_ssymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			float *alpha, int alpha_flag, float *beta,
			int beta_flag, float *a, int lda, float *x, int incx,
			float *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_ssymv{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) float*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  float y_elem;
  float a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_vec;

  float *y_i = y;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *a_i = a;
  float *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }



  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_sdot_testgen(n_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      x_vec, a_vec, seed, &y_elem,
		      &head_r_true_elem, &tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    scopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_sdot_testgen(n_i, i, n_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem = xrand(seed);
      alpha_i[0] = y_elem;
    }
    if (beta_flag == 0) {
      y_elem = xrand(seed);
      beta_i[0] = y_elem;
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_sdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dsymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			double *alpha, int alpha_flag, double *beta,
			int beta_flag, double *a, int lda, double *x,
			int incx, double *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymv{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem;
  double a_elem;
  double x_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  double *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  double *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }



  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_ddot_testgen(n_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      x_vec, a_vec, seed, &y_elem,
		      &head_r_true_elem, &tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    dcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_ddot_testgen(n_i, i, n_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      dsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem = xrand(seed);
      alpha_i[0] = y_elem;
    }
    if (beta_flag == 0) {
      y_elem = xrand(seed);
      beta_i[0] = y_elem;
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_ddot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_csymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymv{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  float y_elem[2];
  float a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *x_vec;

  float *y_i = (float *) y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *x_i = (float *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_cdot_testgen(n_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      x_vec, a_vec, seed, y_elem,
		      head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    ccopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_cdot_testgen(n_i, i, n_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			y_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    ccopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  double a_elem[2];
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *x_i = (double *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_testgen(n_i, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      x_vec, a_vec, seed, y_elem,
		      head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    zcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_testgen(n_i, i, n_i - i, norm,
			blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			y_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = xrand(seed);
      x_elem[1] = xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    zcopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1, x_vec, a_vec, seed,
			y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_csymv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymv_s_s{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  float y_elem[2];
  float a_elem;
  float x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *x_vec;

  float *y_i = (float *) y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_cdot_s_s_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    scopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_cdot_s_s_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      ssy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {

    float *aa_vec;
    float *xx_vec;

    aa_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    xx_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = (float) xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      {
	/* promote to complex */
	int r;
	for (r = 0; r < n_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			xx_vec,
			aa_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

    blas_free(aa_vec);
    blas_free(xx_vec);
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_csymv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymv_s_c{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  float y_elem[2];
  float a_elem;
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *x_vec;

  float *y_i = (float *) y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = a;
  float *x_i = (float *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_cdot_c_s_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    ccopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_cdot_c_s_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      ssy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {

    float *aa_vec;


    aa_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    ccopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      {
	/* promote to complex */
	int r;
	for (r = 0; r < n_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			x_vec,
			aa_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

    blas_free(aa_vec);

  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_csymv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_csymv_c_s{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  float y_elem[2];
  float a_elem[2];
  float x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *x_vec;

  float *y_i = (float *) y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *a_i = (float *) a;
  float *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_cdot_s_c_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    scopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_cdot_s_c_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {


    float *xx_vec;


    xx_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_cdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			xx_vec,
			a_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


    blas_free(xx_vec);
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_d_d{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  double a_elem;
  double x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_d_d_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    dcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_d_d_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      dsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {

    double *aa_vec;
    double *xx_vec;

    aa_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    xx_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = (float) xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      {
	/* promote to complex */
	int r;
	for (r = 0; r < n_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			xx_vec,
			aa_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

    blas_free(aa_vec);
    blas_free(xx_vec);
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_d_z{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  double a_elem;
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = a;
  double *x_i = (double *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_z_d_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    zcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_z_d_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      dsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {

    double *aa_vec;


    aa_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && aa_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }


    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    zcopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      {
	/* promote to complex */
	int r;
	for (r = 0; r < n_i; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			x_vec,
			aa_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

    blas_free(aa_vec);

  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_z_d{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  double a_elem[2];
  double x_elem;
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  double *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  double *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_d_z_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    dcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_d_z_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {


    double *xx_vec;


    xx_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && xx_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = xrand(seed);
      y_elem[1] = xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = xrand(seed);
	a_elem[1] = xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);
    {
      /* promote to complex */
      int r;
      for (r = 0; r < n_i; r++) {
	xx_vec[2 * r] = x_vec[r];
	xx_vec[2 * r + 1] = 0.0;
      }
    }

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_zdot_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			beta, 1,
			xx_vec,
			a_vec,
			seed, y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }


    blas_free(xx_vec);
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dsymv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, float *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymv_s_s{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem;
  float a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  float *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }



  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_ddot_s_s_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, &y_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    scopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_ddot_s_s_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    if (beta_flag == 0) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = (float) xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_ddot_s_s_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dsymv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, double *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymv_s_d{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem;
  float a_elem;
  double x_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  double *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  double *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }



  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_ddot_d_s_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, &y_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    ssy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    dcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_ddot_d_s_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      ssy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    if (beta_flag == 0) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = (float) xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      ssy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_ddot_d_s_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dsymv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, int lda, float *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dsymv_d_s{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem;
  double a_elem;
  float x_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  float *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  float *x_i = x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
  }

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }



  incri = 1;


  incxi = incx;
  incx_veci = 1;



  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_ddot_s_d_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, &y_elem,
			  &head_r_true_elem, &tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;

    /* Copy a_vec to first row of A */
    dsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    scopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_ddot_s_d_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      dsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem = (float) xrand(seed);
      alpha_i[0] = y_elem;
    }
    if (beta_flag == 0) {
      y_elem = (float) xrand(seed);
      beta_i[0] = y_elem;
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;


    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem = (float) xrand(seed);
      x_i[xi] = x_elem;
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_ddot_s_d_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_c_c{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  float a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  float *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  float *x_i = (float *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_c_c_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    ccopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_c_c_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    ccopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_zdot_c_c_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_c_z{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  float a_elem[2];
  double x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  float *a_vec;
  double *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  float *a_i = (float *) a;
  double *x_i = (double *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_z_c_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    csy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    zcopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_z_c_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      csy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    zcopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      csy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_zdot_z_c_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_zsymv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zsymv_z_c{_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  double y_elem[2];
  double a_elem[2];
  float x_elem[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  double *a_vec;
  float *x_vec;

  double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *a_i = (double *) a;
  float *x_i = (float *) x;

  n_i = n;

  /*x_vec, a_vec must have stride of 1 */
  inca = 1;
  inca *= 2;

  a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }

  incyi = incy;
  incyi *= 2;
  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
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

  x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (i = 0; i < n_i; i += incx_veci) {
    x_vec[i] = 0.0;
    x_vec[i + 1] = 0.0;
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */

    BLAS_zdot_c_z_testgen(n_i, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  x_vec, a_vec, seed, y_elem,
			  head_r_true_elem, tail_r_true_elem);

    yi = y_starti;
    ri = 0;
    y_i[yi] = y_elem[0];
    y_i[yi + 1] = y_elem[1];
    head_r_true[ri] = head_r_true_elem[0];
    head_r_true[ri + 1] = head_r_true_elem[1];
    tail_r_true[ri] = tail_r_true_elem[0];
    tail_r_true[ri + 1] = tail_r_true_elem[1];

    /* Copy a_vec to first row of A */
    zsy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    ccopy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for (i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri,
	 yi += incyi) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      BLAS_zdot_c_z_testgen(n_i, i, n_i - i, norm,
			    blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      zsy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }

  } else {







    /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      alpha_i[0] = y_elem[0];
      alpha_i[0 + 1] = y_elem[1];
    }
    if (beta_flag == 0) {
      y_elem[0] = (float) xrand(seed);
      y_elem[1] = (float) xrand(seed);
      beta_i[0] = y_elem[0];
      beta_i[0 + 1] = y_elem[1];
    }


    /*set a, x randomly */

    if (order == blas_colmajor) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incxi = incx;
    incxi *= 2;

    if (incxi < 0) {
      x_starti = (-n + 1) * incxi;
    } else {
      x_starti = 0;
    }

    incai *= 2;
    incaij *= 2;

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem[0] = (float) xrand(seed);
	a_elem[1] = (float) xrand(seed);
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      x_elem[0] = (float) xrand(seed);
      x_elem[1] = (float) xrand(seed);
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* now compute appropriate y vector */

    /* get x */
    ccopy_vector(x_i, n_i, incx, x_vec, 1);


    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      zsy_copy_row(order, uplo, n_i, a, lda, a_vec, i);


      BLAS_zdot_c_z_testgen(n_i, n_i, 0, norm, blas_no_conj, alpha, 1,
			    beta, 1, x_vec, a_vec, seed,
			    y_elem, head_r_true_elem, tail_r_true_elem);

      y_i[yi] = y_elem[0];
      y_i[yi + 1] = y_elem[1];
      head_r_true[ri] = head_r_true_elem[0];
      head_r_true[ri + 1] = head_r_true_elem[1];
      tail_r_true[ri] = tail_r_true_elem[0];
      tail_r_true[ri + 1] = tail_r_true_elem[1];
    }



  }

  blas_free(a_vec);
  blas_free(x_vec);
}
