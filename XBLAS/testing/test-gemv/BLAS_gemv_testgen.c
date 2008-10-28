
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


void BLAS_sgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n,
			float *alpha, int alpha_flag, float *A, int lda,
			float *x, float *beta, int beta_flag, float *y,
			int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) float*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  float y_elem;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_sdot_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed, &y_elem,
		      &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    sge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_dgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n,
			double *alpha, int alpha_flag, double *A, int lda,
			double *x, double *beta, int beta_flag, double *y,
			int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) double*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed, &y_elem,
		      &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    dge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_cgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, void *alpha,
			int alpha_flag, void *A, int lda, void *x, void *beta,
			int beta_flag, void *y, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
		      &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    cge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, void *alpha,
			int alpha_flag, void *A, int lda, void *x, void *beta,
			int beta_flag, void *y, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
		      &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_cgemv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, float *A, int lda,
			    float *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) float*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_s_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    sge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_cgemv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, float *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_c_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    sge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_cgemv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    float *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) float*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  float y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_s_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    cge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, double *A, int lda,
			    double *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) double*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_d_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    dge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, double *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_z_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    dge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    double *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) double*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_d_z_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_dgemv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, float *A, int lda,
			    float *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) float*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_s_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, &y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    sge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_dgemv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, float *A, int lda,
			    double *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) double*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_d_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, &y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    sge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_dgemv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, double *A, int lda,
			    float *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) float*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_s_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, &y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;

    /* copy temp to A */
    dge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_c_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    cge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_z_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    cge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
void BLAS_zgemv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * x            (input/output) void*
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
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  double y_elem[2];

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_c_z_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed, y_elem,
			  &r_true_l[i * incy], &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    /* copy temp to A */
    zge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }

  blas_free(temp);
}
