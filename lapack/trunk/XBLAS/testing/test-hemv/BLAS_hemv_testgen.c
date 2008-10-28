
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void BLAS_sskew_testgen_hemv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo,
			     int n, int randomize,
			     float *alpha, float *beta,
			     float *a, int lda, float *x, int incx, float *y,
			     int incy, int *seed, double *head_r_true,
			     double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
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
 *
 * beta    (input) float*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) float*
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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

  int i, j;
  int yi;
  int aij, ai, ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int incaij, incai;
  int inca_vec;
  int n_i;

  float y_elem;
  float a_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_vec;

  float *y_i = y;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *a_i = a;
  float *x_i = x;

  n_i = n;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }

  incri = 1;


  incx_veci = 1;


  if (randomize == 0) {

    /* Fill in skew matrix A */
    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, yi += incyi) {
      /* x_i has already been copied to x_vec */
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);
      /* skew matricies have zeroed diagonals */
      a_vec[i] = 0.0;
      BLAS_sdot_testgen(n_i, i + 1, n_i - i - 1, norm,
			blas_no_conj, alpha_i, 1,
			beta_i, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);


      sskew_commit_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
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




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	if (i != j) {
	  a_i[aij] = a_elem;
	} else {
	  /* skew matricies have zeroed diagonals */
	  a_i[aij] = 0.0;
	}
      }
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      BLAS_sdot_testgen(n_i, n_i, 0, norm,
			blas_no_conj, alpha_i, 1,
			beta_i, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dskew_testgen_hemv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo,
			     int n, int randomize,
			     double *alpha, double *beta,
			     double *a, int lda, double *x, int incx,
			     double *y, int incy, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
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
 *
 * beta    (input) double*
 *
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) double*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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

  int i, j;
  int yi;
  int aij, ai, ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int incaij, incai;
  int inca_vec;
  int n_i;

  double y_elem;
  double a_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  double *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  double *x_i = x;

  n_i = n;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  dcopy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }

  incri = 1;


  incx_veci = 1;


  if (randomize == 0) {

    /* Fill in skew matrix A */
    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, yi += incyi) {
      /* x_i has already been copied to x_vec */
      dskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);
      /* skew matricies have zeroed diagonals */
      a_vec[i] = 0.0;
      BLAS_ddot_testgen(n_i, i + 1, n_i - i - 1, norm,
			blas_no_conj, alpha_i, 1,
			beta_i, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);


      dskew_commit_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
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




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = xrand(seed);
	if (i != j) {
	  a_i[aij] = a_elem;
	} else {
	  /* skew matricies have zeroed diagonals */
	  a_i[aij] = 0.0;
	}
      }
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      BLAS_ddot_testgen(n_i, n_i, 0, norm,
			blas_no_conj, alpha_i, 1,
			beta_i, 1, x_vec, a_vec, seed,
			&y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dskew_testgen_hemv_d_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo,
				 int n, int randomize,
				 double *alpha, double *beta,
				 double *a, int lda, float *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
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
 *
 * beta    (input) double*
 *
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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

  int i, j;
  int yi;
  int aij, ai, ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int incaij, incai;
  int inca_vec;
  int n_i;

  double y_elem;
  double a_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  float *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *a_i = a;
  float *x_i = x;

  n_i = n;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }

  incri = 1;


  incx_veci = 1;


  if (randomize == 0) {

    /* Fill in skew matrix A */
    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, yi += incyi) {
      /* x_i has already been copied to x_vec */
      dskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);
      /* skew matricies have zeroed diagonals */
      a_vec[i] = 0.0;
      BLAS_ddot_s_d_testgen(n_i, i + 1, n_i - i - 1, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);


      dskew_commit_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
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




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	if (i != j) {
	  a_i[aij] = a_elem;
	} else {
	  /* skew matricies have zeroed diagonals */
	  a_i[aij] = 0.0;
	}
      }
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      dskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      BLAS_ddot_s_d_testgen(n_i, n_i, 0, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dskew_testgen_hemv_s_d(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo,
				 int n, int randomize,
				 double *alpha, double *beta,
				 float *a, int lda, double *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
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
 *
 * beta    (input) double*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) double*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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

  int i, j;
  int yi;
  int aij, ai, ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int incaij, incai;
  int inca_vec;
  int n_i;

  double y_elem;
  float a_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  double *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  double *x_i = x;

  n_i = n;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_vec = (double *) blas_malloc(n_i * sizeof(double));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  dcopy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }

  incri = 1;


  incx_veci = 1;


  if (randomize == 0) {

    /* Fill in skew matrix A */
    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, yi += incyi) {
      /* x_i has already been copied to x_vec */
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);
      /* skew matricies have zeroed diagonals */
      a_vec[i] = 0.0;
      BLAS_ddot_d_s_testgen(n_i, i + 1, n_i - i - 1, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);


      sskew_commit_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
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




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	if (i != j) {
	  a_i[aij] = a_elem;
	} else {
	  /* skew matricies have zeroed diagonals */
	  a_i[aij] = 0.0;
	}
      }
    }

    /* now compute appropriate y vector */

    /* get x */
    dcopy_vector(x_i, n_i, incx, x_vec, 1);

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      BLAS_ddot_d_s_testgen(n_i, n_i, 0, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }
  }

  blas_free(a_vec);
  blas_free(x_vec);
}
void BLAS_dskew_testgen_hemv_s_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo,
				 int n, int randomize,
				 double *alpha, double *beta,
				 float *a, int lda, float *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
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
 *
 * beta    (input) double*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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

  int i, j;
  int yi;
  int aij, ai, ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int incaij, incai;
  int inca_vec;
  int n_i;

  double y_elem;
  float a_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *a_i = a;
  float *x_i = x;

  n_i = n;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n_i; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_vec = (float *) blas_malloc(n_i * sizeof(float));
  if (n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;

  if (incyi < 0) {
    y_starti = (-n + 1) * incyi;
  } else {
    y_starti = 0;
  }

  incri = 1;


  incx_veci = 1;


  if (randomize == 0) {

    /* Fill in skew matrix A */
    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, yi += incyi) {
      /* x_i has already been copied to x_vec */
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);
      /* skew matricies have zeroed diagonals */
      a_vec[i] = 0.0;
      BLAS_ddot_s_s_testgen(n_i, i + 1, n_i - i - 1, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);


      sskew_commit_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      /*commits an element to the generated y */
      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
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




    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	a_elem = (float) xrand(seed);
	if (i != j) {
	  a_i[aij] = a_elem;
	} else {
	  /* skew matricies have zeroed diagonals */
	  a_i[aij] = 0.0;
	}
      }
    }

    /* now compute appropriate y vector */

    /* get x */
    scopy_vector(x_i, n_i, incx, x_vec, 1);

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, ri += incri) {
      sskew_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

      BLAS_ddot_s_s_testgen(n_i, n_i, 0, norm,
			    blas_no_conj, alpha_i, 1,
			    beta_i, 1, x_vec, a_vec, seed,
			    &y_elem, &head_r_true_elem, &tail_r_true_elem);

      y_i[yi] = y_elem;
      head_r_true[ri] = head_r_true_elem;
      tail_r_true[ri] = tail_r_true_elem;
    }
  }

  blas_free(a_vec);
  blas_free(x_vec);
}

void BLAS_chemv_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_chemv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_chemv_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    float *y1;
    float *y2;
    float *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_vec;

    float *y_i = (float *) y;
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;

    float *a_i = (float *) a;
    float *x_i = (float *) x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;
    x_starti *= 2;
    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
      x_vec[i + 1] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_ssymv_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_sskew_testgen_hemv
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a complex vector.  Since x is generated
         as a real vector, we need to perform some scaling.

         There are four cases to consider, depending on the values
         of alpha and beta.

         values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	if (ab == 0) {
	  x_i[xi] = 0.0;
	  x_i[xi + 1] = x0[mi];
	} else {
	  x_i[xi] = x0[mi];
	  x_i[xi + 1] = x0[mi];
	}
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i;
	   i++, yi += incyi, mi += incy1) {
	if (ab == 0) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else if (ab == 2) {
	  y_i[yi] = -2.0 * y2[mi];
	  y_i[yi + 1] = 2.0 * y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in the truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {

	head_r_elem1 = head_r1_true[mi];
	tail_r_elem1 = tail_r1_true[mi];

	head_r_elem2 = head_r2_true[mi];
	tail_r_elem2 = tail_r2_true[mi];

	if (ab == 0) {
	  head_r_true[ri] = -head_r_elem2;
	  tail_r_true[ri] = -tail_r_elem2;
	  head_r_true[ri + 1] = head_r_elem1;
	  tail_r_true[ri + 1] = tail_r_elem1;
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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
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
	  head_r_true[ri] = head_r_elem;
	  tail_r_true[ri] = tail_r_elem;

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
	  head_r_true[ri + 1] = head_r_elem;
	  tail_r_true[ri + 1] = tail_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      float a_elem[2];
      float x_elem[2];
      float y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */


      if (alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      if (beta_flag == 0) {
	beta_i[0] = xrand(seed);
	beta_i[1] = xrand(seed);
      }

      /* Fill in matrix A -- Hermitian. */
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = xrand(seed);
	  a_elem[1] = xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = xrand(seed);
	x_elem[1] = xrand(seed);
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }

      ccopy_vector(x_i, n_i, incx, x_vec, 1);


      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	che_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_cdot_testgen(n_i, n_i, 0,
			  norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec,
			  x_vec, seed,
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
}				/* end BLAS_chemv_testgen */
void BLAS_zhemv_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zhemv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_zhemv_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    double *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    double *x_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    double *x_i = (double *) x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;
    x_starti *= 2;
    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
      x_vec[i + 1] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_dsymv_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_dskew_testgen_hemv
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a complex vector.  Since x is generated
         as a real vector, we need to perform some scaling.

         There are four cases to consider, depending on the values
         of alpha and beta.

         values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	if (ab == 0) {
	  x_i[xi] = 0.0;
	  x_i[xi + 1] = x0[mi];
	} else {
	  x_i[xi] = x0[mi];
	  x_i[xi + 1] = x0[mi];
	}
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i;
	   i++, yi += incyi, mi += incy1) {
	if (ab == 0) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else if (ab == 2) {
	  y_i[yi] = -2.0 * y2[mi];
	  y_i[yi + 1] = 2.0 * y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in the truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {

	head_r_elem1 = head_r1_true[mi];
	tail_r_elem1 = tail_r1_true[mi];

	head_r_elem2 = head_r2_true[mi];
	tail_r_elem2 = tail_r2_true[mi];

	if (ab == 0) {
	  head_r_true[ri] = -head_r_elem2;
	  tail_r_true[ri] = -tail_r_elem2;
	  head_r_true[ri + 1] = head_r_elem1;
	  tail_r_true[ri + 1] = tail_r_elem1;
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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
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
	  head_r_true[ri] = head_r_elem;
	  tail_r_true[ri] = tail_r_elem;

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
	  head_r_true[ri + 1] = head_r_elem;
	  tail_r_true[ri + 1] = tail_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      double a_elem[2];
      double x_elem[2];
      double y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */


      if (alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      if (beta_flag == 0) {
	beta_i[0] = xrand(seed);
	beta_i[1] = xrand(seed);
      }

      /* Fill in matrix A -- Hermitian. */
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = xrand(seed);
	  a_elem[1] = xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = xrand(seed);
	x_elem[1] = xrand(seed);
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }

      zcopy_vector(x_i, n_i, incx, x_vec, 1);


      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	zhe_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_zdot_testgen(n_i, n_i, 0,
			  norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec,
			  x_vec, seed,
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
}				/* end BLAS_zhemv_testgen */
void BLAS_zhemv_c_z_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zhemv_c_z{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_zhemv_c_z_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    double *y1;
    double *y2;
    double *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    double *x_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    float *a_i = (float *) a;
    double *x_i = (double *) x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;
    x_starti *= 2;
    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
      x_vec[i + 1] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_dsymv_s_d_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_dskew_testgen_hemv_s_d
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a complex vector.  Since x is generated
         as a real vector, we need to perform some scaling.

         There are four cases to consider, depending on the values
         of alpha and beta.

         values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	if (ab == 0) {
	  x_i[xi] = 0.0;
	  x_i[xi + 1] = x0[mi];
	} else {
	  x_i[xi] = x0[mi];
	  x_i[xi + 1] = x0[mi];
	}
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i;
	   i++, yi += incyi, mi += incy1) {
	if (ab == 0) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else if (ab == 2) {
	  y_i[yi] = -2.0 * y2[mi];
	  y_i[yi + 1] = 2.0 * y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in the truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {

	head_r_elem1 = head_r1_true[mi];
	tail_r_elem1 = tail_r1_true[mi];

	head_r_elem2 = head_r2_true[mi];
	tail_r_elem2 = tail_r2_true[mi];

	if (ab == 0) {
	  head_r_true[ri] = -head_r_elem2;
	  tail_r_true[ri] = -tail_r_elem2;
	  head_r_true[ri + 1] = head_r_elem1;
	  tail_r_true[ri + 1] = tail_r_elem1;
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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
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
	  head_r_true[ri] = head_r_elem;
	  tail_r_true[ri] = tail_r_elem;

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
	  head_r_true[ri + 1] = head_r_elem;
	  tail_r_true[ri + 1] = tail_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      float a_elem[2];
      double x_elem[2];
      double y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */


      if (alpha_flag == 0) {
	alpha_i[0] = (float) xrand(seed);
	alpha_i[1] = (float) xrand(seed);
      }
      if (beta_flag == 0) {
	beta_i[0] = (float) xrand(seed);
	beta_i[1] = (float) xrand(seed);
      }

      /* Fill in matrix A -- Hermitian. */
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }

      zcopy_vector(x_i, n_i, incx, x_vec, 1);


      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	che_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_zdot_c_z_testgen(n_i, n_i, 0,
			      norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec,
			      x_vec, seed,
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
}				/* end BLAS_zhemv_c_z_testgen */
void BLAS_zhemv_z_c_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zhemv_z_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_zhemv_z_c_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    float *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    float *x_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    float *x_i = (float *) x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;
    x_starti *= 2;
    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
      x_vec[i + 1] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_dsymv_d_s_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_dskew_testgen_hemv_d_s
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a complex vector.  Since x is generated
         as a real vector, we need to perform some scaling.

         There are four cases to consider, depending on the values
         of alpha and beta.

         values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	if (ab == 0) {
	  x_i[xi] = 0.0;
	  x_i[xi + 1] = x0[mi];
	} else {
	  x_i[xi] = x0[mi];
	  x_i[xi + 1] = x0[mi];
	}
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i;
	   i++, yi += incyi, mi += incy1) {
	if (ab == 0) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else if (ab == 2) {
	  y_i[yi] = -2.0 * y2[mi];
	  y_i[yi + 1] = 2.0 * y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in the truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {

	head_r_elem1 = head_r1_true[mi];
	tail_r_elem1 = tail_r1_true[mi];

	head_r_elem2 = head_r2_true[mi];
	tail_r_elem2 = tail_r2_true[mi];

	if (ab == 0) {
	  head_r_true[ri] = -head_r_elem2;
	  tail_r_true[ri] = -tail_r_elem2;
	  head_r_true[ri + 1] = head_r_elem1;
	  tail_r_true[ri + 1] = tail_r_elem1;
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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
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
	  head_r_true[ri] = head_r_elem;
	  tail_r_true[ri] = tail_r_elem;

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
	  head_r_true[ri + 1] = head_r_elem;
	  tail_r_true[ri + 1] = tail_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      double a_elem[2];
      float x_elem[2];
      double y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */


      if (alpha_flag == 0) {
	alpha_i[0] = (float) xrand(seed);
	alpha_i[1] = (float) xrand(seed);
      }
      if (beta_flag == 0) {
	beta_i[0] = (float) xrand(seed);
	beta_i[1] = (float) xrand(seed);
      }

      /* Fill in matrix A -- Hermitian. */
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }

      ccopy_vector(x_i, n_i, incx, x_vec, 1);


      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	zhe_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_zdot_z_c_testgen(n_i, n_i, 0,
			      norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec,
			      x_vec, seed,
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
}				/* end BLAS_zhemv_z_c_testgen */
void BLAS_zhemv_c_c_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zhemv_c_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_zhemv_c_c_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    double *y1;
    double *y2;
    float *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    float *a_i = (float *) a;
    float *x_i = (float *) x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;
    x_starti *= 2;
    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
      x_vec[i + 1] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_dsymv_s_s_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_dskew_testgen_hemv_s_s
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a complex vector.  Since x is generated
         as a real vector, we need to perform some scaling.

         There are four cases to consider, depending on the values
         of alpha and beta.

         values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	if (ab == 0) {
	  x_i[xi] = 0.0;
	  x_i[xi + 1] = x0[mi];
	} else {
	  x_i[xi] = x0[mi];
	  x_i[xi + 1] = x0[mi];
	}
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i;
	   i++, yi += incyi, mi += incy1) {
	if (ab == 0) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else if (ab == 2) {
	  y_i[yi] = -2.0 * y2[mi];
	  y_i[yi + 1] = 2.0 * y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in the truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {

	head_r_elem1 = head_r1_true[mi];
	tail_r_elem1 = tail_r1_true[mi];

	head_r_elem2 = head_r2_true[mi];
	tail_r_elem2 = tail_r2_true[mi];

	if (ab == 0) {
	  head_r_true[ri] = -head_r_elem2;
	  tail_r_true[ri] = -tail_r_elem2;
	  head_r_true[ri + 1] = head_r_elem1;
	  tail_r_true[ri + 1] = tail_r_elem1;
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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
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
	  head_r_true[ri] = head_r_elem;
	  tail_r_true[ri] = tail_r_elem;

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
	  head_r_true[ri + 1] = head_r_elem;
	  tail_r_true[ri + 1] = tail_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      float a_elem[2];
      float x_elem[2];
      double y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */


      if (alpha_flag == 0) {
	alpha_i[0] = (float) xrand(seed);
	alpha_i[1] = (float) xrand(seed);
      }
      if (beta_flag == 0) {
	beta_i[0] = (float) xrand(seed);
	beta_i[1] = (float) xrand(seed);
      }

      /* Fill in matrix A -- Hermitian. */
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = (float) xrand(seed);
	  a_elem[1] = (float) xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem[0] = (float) xrand(seed);
	x_elem[1] = (float) xrand(seed);
	x_i[xi] = x_elem[0];
	x_i[xi + 1] = x_elem[1];
      }

      ccopy_vector(x_i, n_i, incx, x_vec, 1);


      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	che_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_zdot_c_c_testgen(n_i, n_i, 0,
			      norm, blas_no_conj,
			      alpha, 1, beta, 1, a_vec,
			      x_vec, seed,
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
}				/* end BLAS_zhemv_c_c_testgen */
void BLAS_zhemv_z_d_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_zhemv_z_d{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_zhemv_z_d_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    double *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    double *x_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    double *x_i = x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;

    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;


    inca_vec = incx_vec = 1;
    inca_vec *= 2;

    a_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (double *) blas_malloc(n_i * sizeof(double));
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (double *) blas_malloc(n_i * n_i * sizeof(double));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_dsymv_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_dskew_testgen_hemv
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a real vector. 

         There are four cases to consider, depending on the 
         values of alpha and beta.

         values                             scaling
         alpha  beta         alpha    A    x    beta    y     R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	x_i[xi] = x0[mi];
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i; i++, yi += incyi,
	   mi += incy1) {
	if (ab == 1 || ab == 2) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {
	if (ab == 0 || ab == 1) {
	  head_r_true[ri] = head_r1_true[mi];
	  tail_r_true[ri] = tail_r1_true[mi];
	  head_r_true[ri + 1] = head_r2_true[mi];
	  tail_r_true[ri + 1] = tail_r2_true[mi];
	} else if (ab == 2) {
	  head_r_true[ri] = -head_r2_true[mi];
	  tail_r_true[ri] = -tail_r2_true[mi];
	  head_r_true[ri + 1] = head_r1_true[mi];
	  tail_r_true[ri + 1] = tail_r1_true[mi];
	} else {
	  head_r_elem1 = head_r1_true[mi];
	  tail_r_elem1 = tail_r1_true[mi];

	  head_r_elem2 = head_r2_true[mi];
	  tail_r_elem2 = tail_r2_true[mi];

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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      double a_elem[2];
      double x_elem;
      double y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */
      double *xx_vec;
      xx_vec = (double *) blas_malloc(n_i * sizeof(double) * 2);
      if (n_i > 0 && xx_vec == NULL) {
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
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = xrand(seed);
	  a_elem[1] = xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem = xrand(seed);
	x_i[xi] = x_elem;
      }

      dcopy_vector(x_i, n_i, incx, x_vec, 1);

      /* copy the real x_vec into complex xx_vec, so that 
         pure complex test case generator can be called. */
      {
	int k;
	for (k = 0; k < n_i; k++) {
	  xx_vec[2 * k] = x_vec[k];
	  xx_vec[2 * k + 1] = 0.0;
	}
      }

      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	zhe_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_zdot_testgen(n_i, n_i, 0,
			  norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec,
			  xx_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

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
}				/* end BLAS_zhemv_z_d_testgen */
void BLAS_chemv_c_s_testgen(int norm, enum blas_order_type order,
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
 *   Generates the test inputs to BLAS_chemv_c_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hemvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to HEMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
  char *routine_name = "BLAS_chemv_c_s_testgen";
  {

    /* Strategy:  
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incyi, x_starti, y_starti;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int n_i;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    float *y1;
    float *y2;
    float *x0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_vec;

    float *y_i = (float *) y;
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;

    float *a_i = (float *) a;
    float *x_i = x;

    n_i = n;

    ld = n_i;



    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n_i;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n_i;
      inca1ij = incaij = 1;
    }
    if ((0 == incx) || (0 == incy)) {
      BLAS_error(routine_name, 0, 0, NULL);
    }
    if (incx > 0) {
      x_starti = 0;
    } else {
      x_starti = (-n_i + 1) * incx;
    }
    if (incy > 0) {
      y_starti = 0;
    } else {
      y_starti = (-n_i + 1) * incy;
    }
    incri = 1;
    incri *= 2;

    y_starti *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;


    inca_vec = incx_vec = 1;
    inca_vec *= 2;

    a_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
    if (n_i > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * inca_vec; i += inca_vec) {
      a_vec[i] = 0.0;
      a_vec[i + 1] = 0.0;
    }
    x_vec = (float *) blas_malloc(n_i * sizeof(float));
    if (n_i > 0 && x_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    for (i = 0; i < n_i * incx_vec; i += incx_vec) {
      x_vec[i] = 0.0;
    }

    if (randomize == 0) {
      int incx0, incy1, incy2, incmi = 1;
      a1 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      a2 = (float *) blas_malloc(n_i * n_i * sizeof(float));
      if (n_i * n_i > 0 && a2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      for (i = 0; i < n_i * n_i; ++i) {
	a1[i] = a2[i] = 0.0;
      }
      y1 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && y1 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy1 = 1;

      y2 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && y2 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incy2 = 1;

      x0 = (float *) blas_malloc(n_i * sizeof(float));
      if (n_i > 0 && x0 == NULL) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      }
      incx0 = 1;

      for (i = 0; i < n_i; ++i) {
	y1[i] = y2[i] = x0[i] = 0.0;
      }
      head_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r1_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };
      head_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      tail_r2_true = (double *) blas_malloc(n_i * sizeof(double));
      if (n_i > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
	BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
      };

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      /* x0 is output from this call */
      BLAS_ssymv_testgen
	(norm, order, uplo, n, 0, alpha_i, alpha_flag, beta_i,
	 beta_flag, a1, n_i, x0, incx0, y1, incy1, seed, head_r1_true,
	 tail_r1_true);

      /* x0 is now fixed and is input to this call */
      BLAS_sskew_testgen_hemv
	(norm, order, uplo, n, 0, alpha_i, beta_i, a2, n_i, x0,
	 incx0, y2, incy2, seed, head_r2_true, tail_r2_true);


      /* The case where x is a real vector. 

         There are four cases to consider, depending on the 
         values of alpha and beta.

         values                             scaling
         alpha  beta         alpha    A    x    beta    y     R (truth)
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
      for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, a1i += inca1i) {
	for (j = 0, aij = ai, a1ij = a1i; j < n_i;
	     j++, aij += incaij, a1ij += inca1ij) {
	  a_i[aij] = a1[a1ij];
	  a_i[aij + 1] = a2[a1ij];
	}
      }

      /* Fill in x */
      for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi,
	   mi += incx0) {
	x_i[xi] = x0[mi];
      }

      /* Fill in y */
      for (i = 0, yi = y_starti, mi = 0; i < n_i; i++, yi += incyi,
	   mi += incy1) {
	if (ab == 1 || ab == 2) {
	  y_i[yi] = -y2[mi];
	  y_i[yi + 1] = y1[mi];
	} else {
	  y_i[yi] = y1[mi];
	  y_i[yi + 1] = y2[mi];
	}
      }

      /* Fill in truth */
      for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, mi += incmi) {
	if (ab == 0 || ab == 1) {
	  head_r_true[ri] = head_r1_true[mi];
	  tail_r_true[ri] = tail_r1_true[mi];
	  head_r_true[ri + 1] = head_r2_true[mi];
	  tail_r_true[ri + 1] = tail_r2_true[mi];
	} else if (ab == 2) {
	  head_r_true[ri] = -head_r2_true[mi];
	  tail_r_true[ri] = -tail_r2_true[mi];
	  head_r_true[ri + 1] = head_r1_true[mi];
	  tail_r_true[ri + 1] = tail_r1_true[mi];
	} else {
	  head_r_elem1 = head_r1_true[mi];
	  tail_r_elem1 = tail_r1_true[mi];

	  head_r_elem2 = head_r2_true[mi];
	  tail_r_elem2 = tail_r2_true[mi];

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
	  tail_r_true[ri + 1] = tail_r_elem;
	  head_r_true[ri + 1] = head_r_elem;

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
	  tail_r_true[ri] = tail_r_elem;
	  head_r_true[ri] = head_r_elem;
	}
      }


      blas_free(a1);
      blas_free(a2);
      blas_free(y1);
      blas_free(y2);
      blas_free(x0);
      blas_free(head_r1_true);
      blas_free(tail_r1_true);
      blas_free(head_r2_true);
      blas_free(tail_r2_true);
    } else {
      /* get random A, x, then compute y for some cancellation. */
      float a_elem[2];
      float x_elem;
      float y_elem[2];
      double head_r_true_elem[2], tail_r_true_elem[2];

      /* Since mixed real/complex test generator for dot
         scales the vectors, we need to used the non-mixed
         version if x is real (since A is always complex). */
      float *xx_vec;
      xx_vec = (float *) blas_malloc(n_i * sizeof(float) * 2);
      if (n_i > 0 && xx_vec == NULL) {
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
      for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
	for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
	  a_elem[0] = xrand(seed);
	  a_elem[1] = xrand(seed);
	  a_i[aij] = a_elem[0];
	  a_i[aij + 1] = a_elem[1];
	  if (i == j)
	    a_i[aij + 1] = 0.0;
	}
      }

      /* Fill in vector x */
      for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
	x_elem = xrand(seed);
	x_i[xi] = x_elem;
      }

      scopy_vector(x_i, n_i, incx, x_vec, 1);

      /* copy the real x_vec into complex xx_vec, so that 
         pure complex test case generator can be called. */
      {
	int k;
	for (k = 0; k < n_i; k++) {
	  xx_vec[2 * k] = x_vec[k];
	  xx_vec[2 * k + 1] = 0.0;
	}
      }

      for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi,
	   ri += incri) {
	che_copy_row(order, uplo, blas_left_side, n_i, a, lda, a_vec, i);

	BLAS_cdot_testgen(n_i, n_i, 0,
			  norm, blas_no_conj,
			  alpha, 1, beta, 1, a_vec,
			  xx_vec, seed,
			  y_elem, head_r_true_elem, tail_r_true_elem);

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
}				/* end BLAS_chemv_c_s_testgen */
