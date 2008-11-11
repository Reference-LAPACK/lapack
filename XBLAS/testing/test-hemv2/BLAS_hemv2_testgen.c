
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void BLAS_sskmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, float *alpha,
			 float *beta, float *a, int lda, float *x_head,
			 float *x_tail, float *y, int *seed,
			 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) float*
 *
 * beta    (input) float*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) float*
 * x_tail  (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) float*
 *         generated vector y.
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

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  float y_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_head_vec;
  float *x_tail_vec;

  float *y_i = y;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_head_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_head_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_tail_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_tail_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_head_i, n, incx, x_head_vec, 1);
  scopy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;

  yi0 = (incy > 0) ? 0 : -(n - 1) * incyi;

  incri = 1;


  incx_veci = 1;


  /* Fill in skew matrix A */
  for (i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    sskew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    a_vec[i] = 0.0;
    BLAS_sdot2_testgen(n, i + 1, n - i - 1, norm, blas_no_conj, alpha_i, 1,
		       beta_i, 1, x_head_vec, x_tail_vec, a_vec, seed,
		       &y_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;
  }

  blas_free(a_vec);
  blas_free(x_head_vec);
  blas_free(x_tail_vec);
}
void BLAS_dskmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, double *alpha,
			 double *beta, double *a, int lda, double *x_head,
			 double *x_tail, double *y, int *seed,
			 double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) double*
 *
 * beta    (input) double*
 *
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) double*
 * x_tail  (input) double*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) double*
 *         generated vector y.
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

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  double y_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  double *x_head_vec;
  double *x_tail_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_head_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_head_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_tail_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_tail_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  dcopy_vector(x_head_i, n, incx, x_head_vec, 1);
  dcopy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;

  yi0 = (incy > 0) ? 0 : -(n - 1) * incyi;

  incri = 1;


  incx_veci = 1;


  /* Fill in skew matrix A */
  for (i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    dskew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    a_vec[i] = 0.0;
    BLAS_ddot2_testgen(n, i + 1, n - i - 1, norm, blas_no_conj, alpha_i, 1,
		       beta_i, 1, x_head_vec, x_tail_vec, a_vec, seed,
		       &y_elem, &head_r_true_elem, &tail_r_true_elem);

    dskew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;
  }

  blas_free(a_vec);
  blas_free(x_head_vec);
  blas_free(x_tail_vec);
}
void BLAS_dskmv2_testgen_d_s(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, double *a, int lda, float *x_head,
			     float *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) double*
 *
 * beta    (input) double*
 *
 * a       (input/output) double*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) float*
 * x_tail  (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) double*
 *         generated vector y.
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

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  double y_elem;
  double head_r_true_elem, tail_r_true_elem;

  double *a_vec;
  float *x_head_vec;
  float *x_tail_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_head_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_head_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_tail_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_tail_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_head_i, n, incx, x_head_vec, 1);
  scopy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;

  yi0 = (incy > 0) ? 0 : -(n - 1) * incyi;

  incri = 1;


  incx_veci = 1;


  /* Fill in skew matrix A */
  for (i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    dskew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    a_vec[i] = 0.0;
    BLAS_ddot2_s_d_testgen(n, i + 1, n - i - 1, norm, blas_no_conj, alpha_i,
			   1, beta_i, 1, x_head_vec, x_tail_vec, a_vec, seed,
			   &y_elem, &head_r_true_elem, &tail_r_true_elem);

    dskew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;
  }

  blas_free(a_vec);
  blas_free(x_head_vec);
  blas_free(x_tail_vec);
}
void BLAS_dskmv2_testgen_s_d(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, float *a, int lda, double *x_head,
			     double *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) double*
 *
 * beta    (input) double*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) double*
 * x_tail  (input) double*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) double*
 *         generated vector y.
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

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  double y_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  double *x_head_vec;
  double *x_tail_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_head_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_head_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_tail_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_tail_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  dcopy_vector(x_head_i, n, incx, x_head_vec, 1);
  dcopy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;

  yi0 = (incy > 0) ? 0 : -(n - 1) * incyi;

  incri = 1;


  incx_veci = 1;


  /* Fill in skew matrix A */
  for (i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    sskew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    a_vec[i] = 0.0;
    BLAS_ddot2_d_s_testgen(n, i + 1, n - i - 1, norm, blas_no_conj, alpha_i,
			   1, beta_i, 1, x_head_vec, x_tail_vec, a_vec, seed,
			   &y_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;
  }

  blas_free(a_vec);
  blas_free(x_head_vec);
  blas_free(x_tail_vec);
}
void BLAS_dskmv2_testgen_s_s(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, float *a, int lda, float *x_head,
			     float *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) double*
 *
 * beta    (input) double*
 *
 * a       (input/output) float*
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) float*
 * x_tail  (input) float*
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) double*
 *         generated vector y.
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

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  double y_elem;
  double head_r_true_elem, tail_r_true_elem;

  float *a_vec;
  float *x_head_vec;
  float *x_tail_vec;

  double *y_i = y;
  double *alpha_i = alpha;
  double *beta_i = beta;
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;

  /*a_vec must have stride of 1 */
  inca_vec = 1;


  a_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n; i += inca_vec) {
    a_vec[i] = 0.0;
  }
  x_head_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_head_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_tail_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_tail_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  scopy_vector(x_head_i, n, incx, x_head_vec, 1);
  scopy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;

  yi0 = (incy > 0) ? 0 : -(n - 1) * incyi;

  incri = 1;


  incx_veci = 1;


  /* Fill in skew matrix A */
  for (i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    sskew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    a_vec[i] = 0.0;
    BLAS_ddot2_s_s_testgen(n, i + 1, n - i - 1, norm, blas_no_conj, alpha_i,
			   1, beta_i, 1, x_head_vec, x_tail_vec, a_vec, seed,
			   &y_elem, &head_r_true_elem, &tail_r_true_elem);

    sskew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    y_i[yi] = y_elem;
    head_r_true[ri] = head_r_true_elem;
    tail_r_true[ri] = tail_r_true_elem;
  }

  blas_free(a_vec);
  blas_free(x_head_vec);
  blas_free(x_tail_vec);
}


void BLAS_chemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *beta, int beta_flag, void *a,
			 int lda, void *x_head, void *x_tail, void *y,
			 int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_chemv2{_x}
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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) void*
 * x_tail  (input/output) void*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    float *y1;
    float *y2;
    float *x_head_0;
    float *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_head_vec;
    float *x_tail_vec;

    float *y_i = (float *) y;
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;

    float *a_i = (float *) a;
    float *x_head_i = (float *) x_head;
    float *x_tail_i = (float *) x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;
    xi0 *= 2;
    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_ssymv2_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_sskmv2_testgen(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			x_head_0, x_tail_0, y2, seed, head_r2_true,
			tail_r2_true);


    /* The case where x is a complex vector.  Since x is generated
       as a real vector, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    x       beta    y    r (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale r by 1+i, since they are
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      if (ab == 0) {
	x_head_i[xi] = 0.0;
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = 0.0;
	x_tail_i[xi + 1] = x_tail_0[mi];
      } else {
	x_head_i[xi] = x_head_0[mi];
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = x_tail_0[mi];
	x_tail_i[xi + 1] = x_tail_0[mi];
      }
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
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
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_zhemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *beta, int beta_flag, void *a,
			 int lda, void *x_head, void *x_tail, void *y,
			 int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhemv2{_x}
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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) void*
 * x_tail  (input/output) void*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    double *x_head_0;
    double *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    double *x_head_vec;
    double *x_tail_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    double *x_head_i = (double *) x_head;
    double *x_tail_i = (double *) x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;
    xi0 *= 2;
    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_dsymv2_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_dskmv2_testgen(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			x_head_0, x_tail_0, y2, seed, head_r2_true,
			tail_r2_true);


    /* The case where x is a complex vector.  Since x is generated
       as a real vector, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    x       beta    y    r (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale r by 1+i, since they are
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      if (ab == 0) {
	x_head_i[xi] = 0.0;
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = 0.0;
	x_tail_i[xi + 1] = x_tail_0[mi];
      } else {
	x_head_i[xi] = x_head_0[mi];
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = x_tail_0[mi];
	x_tail_i[xi + 1] = x_tail_0[mi];
      }
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
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
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_zhemv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true)

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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) void*
 * x_tail  (input/output) void*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    double *y1;
    double *y2;
    double *x_head_0;
    double *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    double *x_head_vec;
    double *x_tail_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    float *a_i = (float *) a;
    double *x_head_i = (double *) x_head;
    double *x_tail_i = (double *) x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;
    xi0 *= 2;
    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_dsymv2_s_d_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			    x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			    head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_dskmv2_testgen_s_d(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			    x_head_0, x_tail_0, y2, seed, head_r2_true,
			    tail_r2_true);


    /* The case where x is a complex vector.  Since x is generated
       as a real vector, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    x       beta    y    r (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale r by 1+i, since they are
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      if (ab == 0) {
	x_head_i[xi] = 0.0;
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = 0.0;
	x_tail_i[xi + 1] = x_tail_0[mi];
      } else {
	x_head_i[xi] = x_head_0[mi];
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = x_tail_0[mi];
	x_tail_i[xi + 1] = x_tail_0[mi];
      }
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
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
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_zhemv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true)

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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) void*
 * x_tail  (input/output) void*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    float *x_head_0;
    float *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    float *x_head_vec;
    float *x_tail_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    float *x_head_i = (float *) x_head;
    float *x_tail_i = (float *) x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;
    xi0 *= 2;
    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_dsymv2_d_s_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			    x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			    head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_dskmv2_testgen_d_s(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			    x_head_0, x_tail_0, y2, seed, head_r2_true,
			    tail_r2_true);


    /* The case where x is a complex vector.  Since x is generated
       as a real vector, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    x       beta    y    r (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale r by 1+i, since they are
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      if (ab == 0) {
	x_head_i[xi] = 0.0;
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = 0.0;
	x_tail_i[xi + 1] = x_tail_0[mi];
      } else {
	x_head_i[xi] = x_head_0[mi];
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = x_tail_0[mi];
	x_tail_i[xi + 1] = x_tail_0[mi];
      }
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
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
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_zhemv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true)

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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) void*
 * x_tail  (input/output) void*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    double *y1;
    double *y2;
    float *x_head_0;
    float *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_head_vec;
    float *x_tail_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    float *a_i = (float *) a;
    float *x_head_i = (float *) x_head;
    float *x_tail_i = (float *) x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;
    xi0 *= 2;
    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;
    incxi *= 2;

    inca_vec = incx_vec = 1;
    inca_vec *= 2;
    incx_vec *= 2;
    a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_dsymv2_s_s_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			    x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			    head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_dskmv2_testgen_s_s(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			    x_head_0, x_tail_0, y2, seed, head_r2_true,
			    tail_r2_true);


    /* The case where x is a complex vector.  Since x is generated
       as a real vector, we need to perform some scaling.

       There are four cases to consider, depending on the values
       of alpha and beta.

       values                         scaling
       alpha   beta      alpha  A    x       beta    y    r (truth)
       0    1      1                    i               i    i
       1    1      ?                   1+i      1+i         1+i
       2    ?      1         1+i       1+i             2i    2i
       3    ?      ?         1+i       1+i      2i           2i

       Note that we can afford to scale r by 1+i, since they are
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      if (ab == 0) {
	x_head_i[xi] = 0.0;
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = 0.0;
	x_tail_i[xi + 1] = x_tail_0[mi];
      } else {
	x_head_i[xi] = x_head_0[mi];
	x_head_i[xi + 1] = x_head_0[mi];
	x_tail_i[xi] = x_tail_0[mi];
	x_tail_i[xi + 1] = x_tail_0[mi];
      }
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
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
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_zhemv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, double *x_head, double *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true)

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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) double*
 * x_tail  (input/output) double*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    double *a1;
    double *a2;
    double *y1;
    double *y2;
    double *x_head_0;
    double *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    double *a_vec;
    double *x_head_vec;
    double *x_tail_vec;

    double *y_i = (double *) y;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    double *a_i = (double *) a;
    double *x_head_i = x_head;
    double *x_tail_i = x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;

    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;


    inca_vec = incx_vec = 1;
    inca_vec *= 2;

    a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_dsymv2_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_dskmv2_testgen(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			x_head_0, x_tail_0, y2, seed, head_r2_true,
			tail_r2_true);


    /* The case where x is a real vector. 

       There are four cases to consider, depending on the 
       values of alpha and beta.

       values                             scaling
       alpha  beta         alpha    A    x    beta    y     r (truth)
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      x_head_i[xi] = x_head_0[mi];
      x_tail_i[xi] = x_tail_0[mi];
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
      if (ab == 1 || ab == 2) {
	y_i[yi] = -y2[mi];
	y_i[yi + 1] = y1[mi];
      } else {
	y_i[yi] = y1[mi];
	y_i[yi + 1] = y2[mi];
      }
    }

    /* Fill in truth */
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {
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

	/* Set the imaginary part to  r1 + r2 */
	tail_r_true[ri + 1] = tail_r_elem;
	head_r_true[ri + 1] = head_r_elem;

	/* Set the real part to r1 - r2. */
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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}

void BLAS_chemv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, float *x_head, float *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true)

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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * x_head  (input/output) float*
 * x_tail  (input/output) float*
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HEMV2.
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
  {

    /* Strategy:  
       r1 = alpha * A1 * x + beta * y1
       r2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
       and A2 is skew-symmetric.
     */

    int i, j;
    int yi;
    int aij, ai;
    int a1ij, a1i;
    int xi;
    int mi;
    int incx = 1, incy = 1;
    int incyi, xi0, yi0;
    int incaij, incai;
    int inca1ij, inca1i;
    int incxi;
    int inca_vec, incx_vec;
    int ld;
    int ab;
    int ri, incri;

    float *a1;
    float *a2;
    float *y1;
    float *y2;
    float *x_head_0;
    float *x_tail_0;

    double *head_r1_true, *tail_r1_true;
    double *head_r2_true, *tail_r2_true;

    double head_r_elem1, tail_r_elem1;
    double head_r_elem2, tail_r_elem2;
    double head_r_elem, tail_r_elem;

    float *a_vec;
    float *x_head_vec;
    float *x_tail_vec;

    float *y_i = (float *) y;
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;

    float *a_i = (float *) a;
    float *x_head_i = x_head;
    float *x_tail_i = x_tail;

    ld = n;
    if (order == blas_colmajor) {
      inca1i = incai = 1;
      incyi = incy;
      incaij = lda;
      incxi = incx;
      inca1ij = n;
    } else {
      incyi = incy;
      incai = lda;
      incxi = incx;
      inca1i = n;
      inca1ij = incaij = 1;
    }
    xi0 = (incx > 0) ? 0 : -(n - 1) * incx;
    yi0 = (incy > 0) ? 0 : -(n - 1) * incy;
    incri = 1;
    incri *= 2;

    yi0 *= 2;
    incyi *= 2;
    incai *= 2;
    incaij *= 2;


    inca_vec = incx_vec = 1;
    inca_vec *= 2;

    a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && a_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_head_vec = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_head_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_vec = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_tail_vec == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    int incx0, incy1, incy2, incmi = 1;
    a1 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    a2 = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    y1 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && y1 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy1 = 1;

    y2 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && y2 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incy2 = 1;

    x_head_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_head_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    x_tail_0 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x_tail_0 == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }
    incx0 = 1;

    head_r1_true = (double *) blas_malloc(n * sizeof(double));
    tail_r1_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r1_true == NULL || tail_r1_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };
    head_r2_true = (double *) blas_malloc(n * sizeof(double));
    tail_r2_true = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && (head_r2_true == NULL || tail_r2_true == NULL)) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    };

    /* First generate the real portion of A and x. */
    BLAS_ssymv2_testgen(norm, order, uplo, n, alpha_i, alpha_flag, a1, n,
			x_head_0, x_tail_0, beta_i, beta_flag, y1, seed,
			head_r1_true, tail_r1_true);

    /* x0 is now fixed and is input to this call */
    BLAS_sskmv2_testgen(norm, order, uplo, n, alpha_i, beta_i, a2, n,
			x_head_0, x_tail_0, y2, seed, head_r2_true,
			tail_r2_true);


    /* The case where x is a real vector. 

       There are four cases to consider, depending on the 
       values of alpha and beta.

       values                             scaling
       alpha  beta         alpha    A    x    beta    y     r (truth)
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
    for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
      for (j = 0, aij = ai, a1ij = a1i; j < n;
	   j++, aij += incaij, a1ij += inca1ij) {
	a_i[aij] = a1[a1ij];
	a_i[aij + 1] = a2[a1ij];
      }
    }

    /* Fill in x */
    for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
      x_head_i[xi] = x_head_0[mi];
      x_tail_i[xi] = x_tail_0[mi];
    }

    /* Fill in y */
    for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
      if (ab == 1 || ab == 2) {
	y_i[yi] = -y2[mi];
	y_i[yi + 1] = y1[mi];
      } else {
	y_i[yi] = y1[mi];
	y_i[yi + 1] = y2[mi];
      }
    }

    /* Fill in truth */
    for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {
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

	/* Set the imaginary part to  r1 + r2 */
	tail_r_true[ri + 1] = tail_r_elem;
	head_r_true[ri + 1] = head_r_elem;

	/* Set the real part to r1 - r2. */
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
    blas_free(x_head_0);
    blas_free(x_tail_0);
    blas_free(head_r1_true);
    blas_free(tail_r1_true);
    blas_free(head_r2_true);
    blas_free(tail_r2_true);
    blas_free(a_vec);
    blas_free(x_head_vec);
    blas_free(x_tail_vec);
  }
}
