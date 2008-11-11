#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


void BLAS_sgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, float *alpha, int alpha_flag, float *AB,
			 int lda, float *x_head, float *x_tail, float *beta,
			 int beta_flag, float *y, int *seed, double *r_true_l,
			 double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) float*
 * x_tail       (input/output) float*
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
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;
  float *y_i = y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  float y_elem;

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;




  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_sdot2_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
		       seed, &y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;


    /* copy a_vec to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_dgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, double *alpha, int alpha_flag, double *AB,
			 int lda, double *x_head, double *x_tail,
			 double *beta, int beta_flag, double *y, int *seed,
			 double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) double*
 * x_tail       (input/output) double*
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
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem;

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;




  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
		       seed, &y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;


    /* copy a_vec to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_cgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, void *alpha, int alpha_flag, void *AB,
			 int lda, void *x_head, void *x_tail, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  float *x_head_i = (float *) x_head;
  float *x_tail_i = (float *) x_tail;
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  float y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;
  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
		       seed, y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, void *alpha, int alpha_flag, void *AB,
			 int lda, void *x_head, void *x_tail, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  double *x_head_i = (double *) x_head;
  double *x_tail_i = (double *) x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;
  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		       alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
		       seed, y_elem, &r_true_l[i * incy],
		       &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_cgbmv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, float *AB,
			     int lda, float *x_head, float *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) float*
 * x_tail       (input/output) float*
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
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  float y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_s_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];


    /* copy a_vec to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_cgbmv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, float *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  float *x_head_i = (float *) x_head;
  float *x_tail_i = (float *) x_tail;
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  float y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;

  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_c_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];


    /* copy a_vec to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_cgbmv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, float *x_head, float *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) float*
 * x_tail       (input/output) float*
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
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;
  float *y_i = (float *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  float y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_s_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, double *AB,
			     int lda, double *x_head, double *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) double*
 * x_tail       (input/output) double*
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
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_d_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];


    /* copy a_vec to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, double *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  double *x_head_i = (double *) x_head;
  double *x_tail_i = (double *) x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;

  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_z_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];


    /* copy a_vec to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, double *x_head, double *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) double*
 * x_tail       (input/output) double*
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
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_d_z_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_dgbmv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag, float *AB,
			     int lda, float *x_head, float *x_tail,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) float*
 * x_tail       (input/output) float*
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
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem;

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;




  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_s_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;


    /* copy a_vec to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_dgbmv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag, float *AB,
			     int lda, double *x_head, double *x_tail,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) double*
 * x_tail       (input/output) double*
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
  double *x_head_i = x_head;
  double *x_tail_i = x_tail;
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem;

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;




  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_d_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;


    /* copy a_vec to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_dgbmv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag,
			     double *AB, int lda, float *x_head,
			     float *x_tail, double *beta, int beta_flag,
			     double *y, int *seed, double *r_true_l,
			     double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) float*
 * x_tail       (input/output) float*
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
  float *x_head_i = x_head;
  float *x_tail_i = x_tail;
  double *y_i = y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem;

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;




  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_s_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, &y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem;


    /* copy a_vec to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_tail_i[i * incx] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  float *x_head_i = (float *) x_head;
  float *x_tail_i = (float *) x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;
  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_c_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  double *x_head_i = (double *) x_head;
  double *x_tail_i = (double *) x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;
  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_z_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
void BLAS_zgbmv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x_head       (input/output) void*
 * x_tail       (input/output) void*
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
  float *x_head_i = (float *) x_head;
  float *x_tail_i = (float *) x_tail;
  double *y_i = (double *) y;
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *a_vec;
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  double y_elem[2];

  max_mn = MAX(m, n);
  incx = incy = incAB = 1;
  incy *= 2;
  incAB *= 2;
  incx *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to a_vec */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, a_vec, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_c_z_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			   alpha_flag, beta, beta_flag, x_head, x_tail, a_vec,
			   seed, y_elem, &r_true_l[i * incy],
			   &r_true_t[i * incy]);
    y_i[i * incy] = y_elem[0];
    y_i[i * incy + 1] = y_elem[1];

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	a_vec[j + 1] = -a_vec[j + 1];
      }
    }
    /* copy a_vec to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, a_vec, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    x_head_i[i * incx] = 0.0;
    x_head_i[i * incx + 1] = 0.0;
    x_tail_i[i * incx] = 0.0;
    x_tail_i[i * incx + 1] = 0.0;
  }

  blas_free(a_vec);
}
