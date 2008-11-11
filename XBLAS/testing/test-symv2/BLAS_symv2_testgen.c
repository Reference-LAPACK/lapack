
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void BLAS_ssymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, float *alpha,
			 int alpha_flag, float *A, int lda, float *head_x,
			 float *tail_x, float *beta, int beta_flag, float *y,
			 int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
 *
 * beta         (input/output) float*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) float*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *y_i = y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  float y_elem;

  incy = incA = 1;


  temp = (float *) blas_malloc(n * incA * sizeof(float));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    ssy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_sdot2_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, head_x, tail_x, temp,
		       seed, &y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    ssy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_dsymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, double *alpha,
			 int alpha_flag, double *A, int lda, double *head_x,
			 double *tail_x, double *beta, int beta_flag,
			 double *y, int *seed, double *r_true_l,
			 double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem;

  incy = incA = 1;


  temp = (double *) blas_malloc(n * incA * sizeof(double));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    dsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_ddot2_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, head_x, tail_x, temp,
		       seed, &y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    dsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_csymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *A, int lda, void *head_x,
			 void *tail_x, void *beta, int beta_flag, void *y,
			 int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (float *) blas_malloc(n * incA * sizeof(float) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    csy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_cdot2_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, head_x, tail_x, temp,
		       seed, y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    csy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *A, int lda, void *head_x,
			 void *tail_x, void *beta, int beta_flag, void *y,
			 int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (double *) blas_malloc(n * incA * sizeof(double) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    zsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, head_x, tail_x, temp,
		       seed, y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_csymv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, float *A, int lda, float *head_x,
			     float *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;

  temp = (float *) blas_malloc(n * incA * sizeof(float));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    ssy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_cdot2_s_s_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    ssy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_csymv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, float *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;

  temp = (float *) blas_malloc(n * incA * sizeof(float));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    ssy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_cdot2_c_s_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    ssy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_csymv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, float *head_x,
			     float *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (float *) blas_malloc(n * incA * sizeof(float) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    csy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_cdot2_s_c_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    csy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, double *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;

  temp = (double *) blas_malloc(n * incA * sizeof(double));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    dsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_d_d_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    dsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, double *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;

  temp = (double *) blas_malloc(n * incA * sizeof(double));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    dsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_z_d_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    dsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, double *head_x,
			     double *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (double *) blas_malloc(n * incA * sizeof(double) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    zsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_d_z_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_dsymv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, float *A, int lda, float *head_x,
			     float *tail_x, double *beta, int beta_flag,
			     double *y, int *seed, double *r_true_l,
			     double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  double y_elem;

  incy = incA = 1;


  temp = (float *) blas_malloc(n * incA * sizeof(float));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    ssy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_ddot2_s_s_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    ssy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_dsymv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, float *A, int lda,
			     double *head_x, double *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  double y_elem;

  incy = incA = 1;


  temp = (float *) blas_malloc(n * incA * sizeof(float));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    ssy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_ddot2_d_s_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    ssy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_dsymv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, double *A, int lda,
			     float *head_x, float *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem;

  incy = incA = 1;


  temp = (double *) blas_malloc(n * incA * sizeof(double));
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    dsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_ddot2_s_d_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    dsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (float *) blas_malloc(n * incA * sizeof(float) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    csy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_c_c_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    csy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  float *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (float *) blas_malloc(n * incA * sizeof(float) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    csy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_z_c_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    csy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
void BLAS_zsymv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int i;
  double *temp;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;
  temp = (double *) blas_malloc(n * incA * sizeof(double) * 2);
  if (n * incA > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1;
      beta_flag = 1;
    }

    zsy_copy_row(order, uplo, n, A, lda, temp, i);
    BLAS_zdot2_c_z_testgen(n, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, head_x, tail_x, temp,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zsy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  blas_free(temp);
}
