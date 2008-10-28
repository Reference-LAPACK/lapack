
#include "blas_extended.h"
#include "blas_extended_test.h"


void BLAS_stpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *tp, float *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) float*
 *
 * x            (input/output) float*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  float *temp;
  float beta_zero_fake;
  int beta_flag = 1;
  float r_zero_fake;
  float one_this;
  float tp_swapTemp1;
  float tp_swapTemp2;
  float x_swapTemp1;
  float x_swapTemp2;
  double head_r_true_swapTemp1, tail_r_true_swapTemp1;
  double head_r_true_swapTemp2, tail_r_true_swapTemp2;

  float *x_i = x;
  float *alpha_i = alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake = 0.0;
  r_zero_fake = 0.0;
  one_this = 1.0;




  temp = (float *) blas_malloc(n * 2 * inc_index * sizeof(float));
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    tp[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    *alpha_i = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_sdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			alpha, beta_flag, x, temp, seed, &r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_sdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			&beta_zero_fake, beta_flag, x, temp, seed,
			&r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    }

/* one type is s */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this;
      x_i[length * inc_index] = r_zero_fake;
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1 = temp[(n - 1 - j) * inctp];
	tp_swapTemp2 = temp[j * inctp];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2;
	temp[j * inctp] = tp_swapTemp1;
      }

      stpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      stpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1 = x_i[(n - 1 - j) * inc_index];
      x_swapTemp2 = x_i[j * inc_index];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2;
      x_i[j * inc_index] = x_swapTemp1;

      head_r_true_swapTemp1 = head_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1 = tail_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp2 = head_r_true[j * inc_index];
      tail_r_true_swapTemp2 = tail_r_true[j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2;
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2;
      head_r_true[j * inc_index] = head_r_true_swapTemp1;
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1;
    }
  }

  blas_free(temp);
}

void BLAS_dtpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *tp, double *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) double*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  double *temp;
  double beta_zero_fake;
  int beta_flag = 1;
  double r_zero_fake;
  double one_this;
  double tp_swapTemp1;
  double tp_swapTemp2;
  double x_swapTemp1;
  double x_swapTemp2;
  double head_r_true_swapTemp1, tail_r_true_swapTemp1;
  double head_r_true_swapTemp2, tail_r_true_swapTemp2;

  double *x_i = x;
  double *alpha_i = alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake = 0.0;
  r_zero_fake = 0.0;
  one_this = 1.0;




  temp = (double *) blas_malloc(n * 2 * inc_index * sizeof(double));
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    tp[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    *alpha_i = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_ddot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			alpha, beta_flag, x, temp, seed, &r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_ddot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			&beta_zero_fake, beta_flag, x, temp, seed,
			&r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    }

/* one type is d */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this;
      x_i[length * inc_index] = r_zero_fake;
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1 = temp[(n - 1 - j) * inctp];
	tp_swapTemp2 = temp[j * inctp];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2;
	temp[j * inctp] = tp_swapTemp1;
      }

      dtpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      dtpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1 = x_i[(n - 1 - j) * inc_index];
      x_swapTemp2 = x_i[j * inc_index];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2;
      x_i[j * inc_index] = x_swapTemp1;

      head_r_true_swapTemp1 = head_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1 = tail_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp2 = head_r_true[j * inc_index];
      tail_r_true_swapTemp2 = tail_r_true[j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2;
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2;
      head_r_true[j * inc_index] = head_r_true_swapTemp1;
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1;
    }
  }

  blas_free(temp);
}

void BLAS_ctpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *tp, void *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  float *temp;
  float beta_zero_fake[2];
  int beta_flag = 1;
  float r_zero_fake[2];
  float one_this[2];
  float tp_swapTemp1[2];
  float tp_swapTemp2[2];
  float x_swapTemp1[2];
  float x_swapTemp2[2];
  double head_r_true_swapTemp1[2], tail_r_true_swapTemp1[2];
  double head_r_true_swapTemp2[2], tail_r_true_swapTemp2[2];

  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;
  r_zero_fake[0] = r_zero_fake[1] = 0.0;
  one_this[0] = 1.0;
  one_this[1] = 0.0;

  inctp *= 2;
  inc_index *= 2;

  temp = (float *) blas_malloc(n * 2 * inc_index * sizeof(float) * 2);
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    ((float *) tp)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
      temp[j + 1] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_cdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			alpha, beta_flag, x, temp, seed, r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_cdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			beta_zero_fake, beta_flag, x, temp, seed, r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    }

/* one type is c */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this[0];
      temp[length * inctp + 1] = one_this[1];
      x_i[length * inc_index] = r_zero_fake[0];
      x_i[length * inc_index + 1] = r_zero_fake[1];
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1[0] = temp[(n - 1 - j) * inctp];
	tp_swapTemp1[1] = temp[(n - 1 - j) * inctp + 1];
	tp_swapTemp2[0] = temp[j * inctp];
	tp_swapTemp2[1] = temp[j * inctp + 1];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2[0];
	temp[(n - 1 - j) * inctp + 1] = tp_swapTemp2[1];
	temp[j * inctp] = tp_swapTemp1[0];
	temp[j * inctp + 1] = tp_swapTemp1[1];
      }

      ctpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      ctpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1[0] = x_i[(n - 1 - j) * inc_index];
      x_swapTemp1[1] = x_i[(n - 1 - j) * inc_index + 1];
      x_swapTemp2[0] = x_i[j * inc_index];
      x_swapTemp2[1] = x_i[j * inc_index + 1];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2[0];
      x_i[(n - 1 - j) * inc_index + 1] = x_swapTemp2[1];
      x_i[j * inc_index] = x_swapTemp1[0];
      x_i[j * inc_index + 1] = x_swapTemp1[1];

      head_r_true_swapTemp1[0] = head_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp1[1] = head_r_true[1 + (n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[0] = tail_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[1] = tail_r_true[1 + (n - 1 - j) * inc_index];
      head_r_true_swapTemp2[0] = head_r_true[j * inc_index];
      head_r_true_swapTemp2[1] = head_r_true[1 + j * inc_index];
      tail_r_true_swapTemp2[0] = tail_r_true[j * inc_index];
      tail_r_true_swapTemp2[1] = tail_r_true[1 + j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2[0];
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2[0];
      head_r_true[1 + (n - 1 - j) * inc_index] = head_r_true_swapTemp2[1];
      tail_r_true[1 + (n - 1 - j) * inc_index] = tail_r_true_swapTemp2[1];
      head_r_true[j * inc_index] = head_r_true_swapTemp1[0];
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1[0];
      head_r_true[1 + j * inc_index] = head_r_true_swapTemp1[1];
      tail_r_true[1 + j * inc_index] = tail_r_true_swapTemp1[1];
    }
  }

  blas_free(temp);
}

void BLAS_ztpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *tp, void *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  double *temp;
  double beta_zero_fake[2];
  int beta_flag = 1;
  double r_zero_fake[2];
  double one_this[2];
  double tp_swapTemp1[2];
  double tp_swapTemp2[2];
  double x_swapTemp1[2];
  double x_swapTemp2[2];
  double head_r_true_swapTemp1[2], tail_r_true_swapTemp1[2];
  double head_r_true_swapTemp2[2], tail_r_true_swapTemp2[2];

  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;
  r_zero_fake[0] = r_zero_fake[1] = 0.0;
  one_this[0] = 1.0;
  one_this[1] = 0.0;

  inctp *= 2;
  inc_index *= 2;

  temp = (double *) blas_malloc(n * 2 * inc_index * sizeof(double) * 2);
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    ((double *) tp)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
      temp[j + 1] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_zdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			alpha, beta_flag, x, temp, seed, r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_zdot_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			alpha, alpha_flag,
			beta_zero_fake, beta_flag, x, temp, seed, r_zero_fake,
			&((double *) head_r_true)[i * inc_index],
			&((double *) tail_r_true)[i * inc_index]);
    }

/* one type is z */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this[0];
      temp[length * inctp + 1] = one_this[1];
      x_i[length * inc_index] = r_zero_fake[0];
      x_i[length * inc_index + 1] = r_zero_fake[1];
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1[0] = temp[(n - 1 - j) * inctp];
	tp_swapTemp1[1] = temp[(n - 1 - j) * inctp + 1];
	tp_swapTemp2[0] = temp[j * inctp];
	tp_swapTemp2[1] = temp[j * inctp + 1];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2[0];
	temp[(n - 1 - j) * inctp + 1] = tp_swapTemp2[1];
	temp[j * inctp] = tp_swapTemp1[0];
	temp[j * inctp + 1] = tp_swapTemp1[1];
      }

      ztpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      ztpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1[0] = x_i[(n - 1 - j) * inc_index];
      x_swapTemp1[1] = x_i[(n - 1 - j) * inc_index + 1];
      x_swapTemp2[0] = x_i[j * inc_index];
      x_swapTemp2[1] = x_i[j * inc_index + 1];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2[0];
      x_i[(n - 1 - j) * inc_index + 1] = x_swapTemp2[1];
      x_i[j * inc_index] = x_swapTemp1[0];
      x_i[j * inc_index + 1] = x_swapTemp1[1];

      head_r_true_swapTemp1[0] = head_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp1[1] = head_r_true[1 + (n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[0] = tail_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[1] = tail_r_true[1 + (n - 1 - j) * inc_index];
      head_r_true_swapTemp2[0] = head_r_true[j * inc_index];
      head_r_true_swapTemp2[1] = head_r_true[1 + j * inc_index];
      tail_r_true_swapTemp2[0] = tail_r_true[j * inc_index];
      tail_r_true_swapTemp2[1] = tail_r_true[1 + j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2[0];
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2[0];
      head_r_true[1 + (n - 1 - j) * inc_index] = head_r_true_swapTemp2[1];
      tail_r_true[1 + (n - 1 - j) * inc_index] = tail_r_true_swapTemp2[1];
      head_r_true[j * inc_index] = head_r_true_swapTemp1[0];
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1[0];
      head_r_true[1 + j * inc_index] = head_r_true_swapTemp1[1];
      tail_r_true[1 + j * inc_index] = tail_r_true_swapTemp1[1];
    }
  }

  blas_free(temp);
}

void BLAS_dtpmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *tp, double *x, int *seed,
			  double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) float*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  float *temp;
  double beta_zero_fake;
  int beta_flag = 1;
  double r_zero_fake;
  float one_this;
  float tp_swapTemp1;
  float tp_swapTemp2;
  double x_swapTemp1;
  double x_swapTemp2;
  double head_r_true_swapTemp1, tail_r_true_swapTemp1;
  double head_r_true_swapTemp2, tail_r_true_swapTemp2;

  double *x_i = x;
  double *alpha_i = alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake = 0.0;
  r_zero_fake = 0.0;
  one_this = 1.0;




  temp = (float *) blas_malloc(n * 2 * inc_index * sizeof(float));
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    tp[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    *alpha_i = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_ddot_d_s_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    alpha, beta_flag, x, temp, seed, &r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_ddot_d_s_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    &beta_zero_fake, beta_flag, x, temp, seed,
			    &r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    }

/* one type is d */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this;
      x_i[length * inc_index] = r_zero_fake;
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1 = temp[(n - 1 - j) * inctp];
	tp_swapTemp2 = temp[j * inctp];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2;
	temp[j * inctp] = tp_swapTemp1;
      }

      stpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      stpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1 = x_i[(n - 1 - j) * inc_index];
      x_swapTemp2 = x_i[j * inc_index];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2;
      x_i[j * inc_index] = x_swapTemp1;

      head_r_true_swapTemp1 = head_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1 = tail_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp2 = head_r_true[j * inc_index];
      tail_r_true_swapTemp2 = tail_r_true[j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2;
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2;
      head_r_true[j * inc_index] = head_r_true_swapTemp1;
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1;
    }
  }

  blas_free(temp);
}

void BLAS_ztpmv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  float *temp;
  double beta_zero_fake[2];
  int beta_flag = 1;
  double r_zero_fake[2];
  float one_this[2];
  float tp_swapTemp1[2];
  float tp_swapTemp2[2];
  double x_swapTemp1[2];
  double x_swapTemp2[2];
  double head_r_true_swapTemp1[2], tail_r_true_swapTemp1[2];
  double head_r_true_swapTemp2[2], tail_r_true_swapTemp2[2];

  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;
  r_zero_fake[0] = r_zero_fake[1] = 0.0;
  one_this[0] = 1.0;
  one_this[1] = 0.0;

  inctp *= 2;
  inc_index *= 2;

  temp = (float *) blas_malloc(n * 2 * inc_index * sizeof(float) * 2);
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    ((float *) tp)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
      temp[j + 1] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_zdot_z_c_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    alpha, beta_flag, x, temp, seed, r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_zdot_z_c_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    beta_zero_fake, beta_flag, x, temp, seed,
			    r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    }

/* one type is z */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this[0];
      temp[length * inctp + 1] = one_this[1];
      x_i[length * inc_index] = r_zero_fake[0];
      x_i[length * inc_index + 1] = r_zero_fake[1];
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1[0] = temp[(n - 1 - j) * inctp];
	tp_swapTemp1[1] = temp[(n - 1 - j) * inctp + 1];
	tp_swapTemp2[0] = temp[j * inctp];
	tp_swapTemp2[1] = temp[j * inctp + 1];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2[0];
	temp[(n - 1 - j) * inctp + 1] = tp_swapTemp2[1];
	temp[j * inctp] = tp_swapTemp1[0];
	temp[j * inctp + 1] = tp_swapTemp1[1];
      }

      ctpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      ctpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1[0] = x_i[(n - 1 - j) * inc_index];
      x_swapTemp1[1] = x_i[(n - 1 - j) * inc_index + 1];
      x_swapTemp2[0] = x_i[j * inc_index];
      x_swapTemp2[1] = x_i[j * inc_index + 1];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2[0];
      x_i[(n - 1 - j) * inc_index + 1] = x_swapTemp2[1];
      x_i[j * inc_index] = x_swapTemp1[0];
      x_i[j * inc_index + 1] = x_swapTemp1[1];

      head_r_true_swapTemp1[0] = head_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp1[1] = head_r_true[1 + (n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[0] = tail_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[1] = tail_r_true[1 + (n - 1 - j) * inc_index];
      head_r_true_swapTemp2[0] = head_r_true[j * inc_index];
      head_r_true_swapTemp2[1] = head_r_true[1 + j * inc_index];
      tail_r_true_swapTemp2[0] = tail_r_true[j * inc_index];
      tail_r_true_swapTemp2[1] = tail_r_true[1 + j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2[0];
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2[0];
      head_r_true[1 + (n - 1 - j) * inc_index] = head_r_true_swapTemp2[1];
      tail_r_true[1 + (n - 1 - j) * inc_index] = tail_r_true_swapTemp2[1];
      head_r_true[j * inc_index] = head_r_true_swapTemp1[0];
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1[0];
      head_r_true[1 + j * inc_index] = head_r_true_swapTemp1[1];
      tail_r_true[1 + j * inc_index] = tail_r_true_swapTemp1[1];
    }
  }

  blas_free(temp);
}

void BLAS_ctpmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) float*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  float *temp;
  float beta_zero_fake[2];
  int beta_flag = 1;
  float r_zero_fake[2];
  float one_this;
  float tp_swapTemp1;
  float tp_swapTemp2;
  float x_swapTemp1[2];
  float x_swapTemp2[2];
  double head_r_true_swapTemp1[2], tail_r_true_swapTemp1[2];
  double head_r_true_swapTemp2[2], tail_r_true_swapTemp2[2];

  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;
  r_zero_fake[0] = r_zero_fake[1] = 0.0;
  one_this = 1.0;


  inc_index *= 2;

  temp = (float *) blas_malloc(n * 2 * inc_index * sizeof(float));
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    tp[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_cdot_c_s_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    alpha, beta_flag, x, temp, seed, r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_cdot_c_s_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    beta_zero_fake, beta_flag, x, temp, seed,
			    r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    }

/* one type is c */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this;
      x_i[length * inc_index] = r_zero_fake[0];
      x_i[length * inc_index + 1] = r_zero_fake[1];
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1 = temp[(n - 1 - j) * inctp];
	tp_swapTemp2 = temp[j * inctp];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2;
	temp[j * inctp] = tp_swapTemp1;
      }

      stpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      stpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1[0] = x_i[(n - 1 - j) * inc_index];
      x_swapTemp1[1] = x_i[(n - 1 - j) * inc_index + 1];
      x_swapTemp2[0] = x_i[j * inc_index];
      x_swapTemp2[1] = x_i[j * inc_index + 1];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2[0];
      x_i[(n - 1 - j) * inc_index + 1] = x_swapTemp2[1];
      x_i[j * inc_index] = x_swapTemp1[0];
      x_i[j * inc_index + 1] = x_swapTemp1[1];

      head_r_true_swapTemp1[0] = head_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp1[1] = head_r_true[1 + (n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[0] = tail_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[1] = tail_r_true[1 + (n - 1 - j) * inc_index];
      head_r_true_swapTemp2[0] = head_r_true[j * inc_index];
      head_r_true_swapTemp2[1] = head_r_true[1 + j * inc_index];
      tail_r_true_swapTemp2[0] = tail_r_true[j * inc_index];
      tail_r_true_swapTemp2[1] = tail_r_true[1 + j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2[0];
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2[0];
      head_r_true[1 + (n - 1 - j) * inc_index] = head_r_true_swapTemp2[1];
      tail_r_true[1 + (n - 1 - j) * inc_index] = tail_r_true_swapTemp2[1];
      head_r_true[j * inc_index] = head_r_true_swapTemp1[0];
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1[0];
      head_r_true[1 + j * inc_index] = head_r_true_swapTemp1[1];
      tail_r_true[1 + j * inc_index] = tail_r_true_swapTemp1[1];
    }
  }

  blas_free(temp);
}

void BLAS_ztpmv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) double*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;
  double *temp;
  double beta_zero_fake[2];
  int beta_flag = 1;
  double r_zero_fake[2];
  double one_this;
  double tp_swapTemp1;
  double tp_swapTemp2;
  double x_swapTemp1[2];
  double x_swapTemp2[2];
  double head_r_true_swapTemp1[2], tail_r_true_swapTemp1[2];
  double head_r_true_swapTemp2[2], tail_r_true_swapTemp2[2];

  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  int inctp = 1;
  int inc_index = 1;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;
  r_zero_fake[0] = r_zero_fake[1] = 0.0;
  one_this = 1.0;


  inc_index *= 2;

  temp = (double *) blas_malloc(n * 2 * inc_index * sizeof(double));
  if (n * 2 * inc_index > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < n * 2 * inc_index; i += inc_index) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (n - 1 + n - 1 + 1) * n * 2 * inctp; i++) {
    tp[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
    alpha_flag = 1;
  }

  for (i = 0; i < n; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else {
      n_mix = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j = 0; j < n * inc_index; j += inc_index) {
      temp[j] = 0.0;
    }
    length = i + 1;

    if (diag == blas_unit_diag) {
      length = i;


      BLAS_zdot_z_d_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    alpha, beta_flag, x, temp, seed, r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    } else {
      BLAS_zdot_z_d_testgen(length, n_fix2, n_mix, norm, blas_no_conj,
			    alpha, alpha_flag,
			    beta_zero_fake, beta_flag, x, temp, seed,
			    r_zero_fake,
			    &((double *) head_r_true)[i * inc_index],
			    &((double *) tail_r_true)[i * inc_index]);
    }

/* one type is z */
    if (diag == blas_unit_diag) {
      temp[length * inctp] = one_this;
      x_i[length * inc_index] = r_zero_fake[0];
      x_i[length * inc_index + 1] = r_zero_fake[1];
    }


    /* copy temp to tp */
    if (((order == blas_rowmajor) && (trans == blas_no_trans)
	 && (uplo == blas_upper))
	|| ((order == blas_colmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_rowmajor) && (trans != blas_no_trans)
	    && (uplo == blas_lower))
	|| ((order == blas_colmajor) && (trans == blas_no_trans)
	    && (uplo == blas_upper))) {
      for (j = 0; j < n / 2; j++) {
	tp_swapTemp1 = temp[(n - 1 - j) * inctp];
	tp_swapTemp2 = temp[j * inctp];
	temp[(n - 1 - j) * inctp] = tp_swapTemp2;
	temp[j * inctp] = tp_swapTemp1;
      }

      dtpmv_commit_row(order, uplo, trans, n, tp, temp, (n - 1 - i));
    } else
      if (((order == blas_rowmajor) && (trans == blas_no_trans)
	   && (uplo == blas_lower))
	  || ((order == blas_colmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_rowmajor) && (trans != blas_no_trans)
	      && (uplo == blas_upper))
	  || ((order == blas_colmajor) && (trans == blas_no_trans)
	      && (uplo == blas_lower))) {
      dtpmv_commit_row(order, uplo, trans, n, tp, temp, i);
    }

  }

  /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
  if (((order == blas_rowmajor) && (trans == blas_no_trans)
       && (uplo == blas_upper))
      || ((order == blas_colmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_rowmajor) && (trans != blas_no_trans)
	  && (uplo == blas_lower))
      || ((order == blas_colmajor) && (trans == blas_no_trans)
	  && (uplo == blas_upper))) {
    for (j = 0; j < n / 2; j++) {
      x_swapTemp1[0] = x_i[(n - 1 - j) * inc_index];
      x_swapTemp1[1] = x_i[(n - 1 - j) * inc_index + 1];
      x_swapTemp2[0] = x_i[j * inc_index];
      x_swapTemp2[1] = x_i[j * inc_index + 1];
      x_i[(n - 1 - j) * inc_index] = x_swapTemp2[0];
      x_i[(n - 1 - j) * inc_index + 1] = x_swapTemp2[1];
      x_i[j * inc_index] = x_swapTemp1[0];
      x_i[j * inc_index + 1] = x_swapTemp1[1];

      head_r_true_swapTemp1[0] = head_r_true[(n - 1 - j) * inc_index];
      head_r_true_swapTemp1[1] = head_r_true[1 + (n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[0] = tail_r_true[(n - 1 - j) * inc_index];
      tail_r_true_swapTemp1[1] = tail_r_true[1 + (n - 1 - j) * inc_index];
      head_r_true_swapTemp2[0] = head_r_true[j * inc_index];
      head_r_true_swapTemp2[1] = head_r_true[1 + j * inc_index];
      tail_r_true_swapTemp2[0] = tail_r_true[j * inc_index];
      tail_r_true_swapTemp2[1] = tail_r_true[1 + j * inc_index];
      head_r_true[(n - 1 - j) * inc_index] = head_r_true_swapTemp2[0];
      tail_r_true[(n - 1 - j) * inc_index] = tail_r_true_swapTemp2[0];
      head_r_true[1 + (n - 1 - j) * inc_index] = head_r_true_swapTemp2[1];
      tail_r_true[1 + (n - 1 - j) * inc_index] = tail_r_true_swapTemp2[1];
      head_r_true[j * inc_index] = head_r_true_swapTemp1[0];
      tail_r_true[j * inc_index] = tail_r_true_swapTemp1[0];
      head_r_true[1 + j * inc_index] = head_r_true_swapTemp1[1];
      tail_r_true[1 + j * inc_index] = tail_r_true_swapTemp1[1];
    }
  }

  blas_free(temp);
}
