#include <stdio.h>
#include <assert.h>
#include "blas_extended.h"
#include "blas_extended_test.h"


void BLAS_sdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, float *alpha,
			int alpha_flag, float *beta, int beta_flag,
			float *head_x, float *tail_x, float *y, int *seed,
			float *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_sdot2{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) float*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) float*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) float*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float *alpha_i = alpha;
  float *beta_i = beta;
  float *r_i = r;
  float *head_x_i = head_x;
  float *tail_x_i = tail_x;
  float *y_i = y;
  float alpha_tmp;
  float beta_tmp;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp = *alpha_i;
  beta_tmp = *beta_i;

  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_tmp, alpha_flag,
		     &beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     r_true_l, r_true_t);

  *alpha_i = alpha_tmp;
  *beta_i = beta_tmp;
  *r_i = r_tmp;
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_sdot2_testgen */

void BLAS_ddot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, double *alpha,
			int alpha_flag, double *beta, int beta_flag,
			double *head_x, double *tail_x, double *y, int *seed,
			double *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot2{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *r_i = r;
  double *head_x_i = head_x;
  double *tail_x_i = tail_x;
  double *y_i = y;
  double alpha_tmp;
  double beta_tmp;
  double r_tmp;
  double *head_x_vec;
  double *tail_x_vec;
  double *y_vec;

  alpha_tmp = *alpha_i;
  beta_tmp = *beta_i;

  head_x_vec = (double *) blas_malloc(3 * n * sizeof(double));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_tmp, alpha_flag,
		     &beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     r_true_l, r_true_t);

  *alpha_i = alpha_tmp;
  *beta_i = beta_tmp;
  *r_i = r_tmp;
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_ddot2_testgen */

void BLAS_cdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot2{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *r_i = (float *) r;
  float *head_x_i = (float *) head_x;
  float *tail_x_i = (float *) tail_x;
  float *y_i = (float *) y;
  float alpha_tmp[2];
  float beta_tmp[2];
  float r_tmp[2];
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp[0] = alpha_i[0];
  alpha_tmp[1] = alpha_i[1];
  beta_tmp[0] = beta_i[0];
  beta_tmp[1] = beta_i[1];
  inc *= 2;
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float) * 2);
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
    y_vec[i] = y_i[i];
    y_vec[i + 1] = y_i[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot2(n, n_fix2, n_mix, norm, conj,
		     alpha_tmp, alpha_flag,
		     beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, r_tmp,
		     r_true_l, r_true_t);

  alpha_i[0] = alpha_tmp[0];
  alpha_i[1] = alpha_tmp[1];
  beta_i[0] = beta_tmp[0];
  beta_i[1] = beta_tmp[1];
  r_i[0] = r_tmp[0];
  r_i[1] = r_tmp[1];
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    head_x_i[i + 1] = head_x_vec[i + 1];
    tail_x_i[i] = tail_x_vec[i];
    tail_x_i[i + 1] = tail_x_vec[i + 1];
    y_i[i] = y_vec[i];
    y_i[i + 1] = y_vec[i + 1];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_cdot2_testgen */

void BLAS_zdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *head_x_i = (double *) head_x;
  double *tail_x_i = (double *) tail_x;
  double *y_i = (double *) y;
  double alpha_tmp[2];
  double beta_tmp[2];
  double r_tmp[2];
  double *head_x_vec;
  double *tail_x_vec;
  double *y_vec;

  alpha_tmp[0] = alpha_i[0];
  alpha_tmp[1] = alpha_i[1];
  beta_tmp[0] = beta_i[0];
  beta_tmp[1] = beta_i[1];
  inc *= 2;
  head_x_vec = (double *) blas_malloc(3 * n * sizeof(double) * 2);
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
    y_vec[i] = y_i[i];
    y_vec[i + 1] = y_i[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_zdot2(n, n_fix2, n_mix, norm, conj,
		     alpha_tmp, alpha_flag,
		     beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, r_tmp,
		     r_true_l, r_true_t);

  alpha_i[0] = alpha_tmp[0];
  alpha_i[1] = alpha_tmp[1];
  beta_i[0] = beta_tmp[0];
  beta_i[1] = beta_tmp[1];
  r_i[0] = r_tmp[0];
  r_i[1] = r_tmp[1];
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    head_x_i[i + 1] = head_x_vec[i + 1];
    tail_x_i[i] = tail_x_vec[i];
    tail_x_i[i + 1] = tail_x_vec[i + 1];
    y_i[i] = y_vec[i];
    y_i[i + 1] = y_vec[i + 1];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_testgen */

void BLAS_cdot2_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    float *head_x, float *tail_x, float *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot2_s_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *r_i = (float *) r;
  float *head_x_i = head_x;
  float *tail_x_i = tail_x;
  float *y_i = y;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 && ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.))) {	/* beta_i == 0 or 1 */
	r_i[0] = r_tmp;
	r_i[1] = 0.0;
      } else {			/* beta_i *= (1-i), r_i *= (1+i)/2 --> prod = 1 */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = r_tmp / 2.;
	r_i[1] = r_tmp / 2.;
      }
      r_true_l[1] = r_true_t[1] = 0.0;
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i *= (1+i) */
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    } else {			/* beta_i *= (1-i), r_i *= i */
      beta_i[0] = beta_i_r;
      beta_i[1] = -beta_i_r;
      r_i[0] = 0.;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = -beta_i_r;
    r_i[0] = 0;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  }
  for (i = 0; i < n; ++i) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }


  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_cdot2_s_s_testgen */

void BLAS_cdot2_s_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    float *head_x, float *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot2_s_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *r_i = (float *) r;
  float *head_x_i = head_x;
  float *tail_x_i = tail_x;
  float *y_i = (float *) y;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[2 * i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i is 0 or 1 --> r_i *= 2i --> prod = 2i */
	r_i[0] = 0.0;
	r_i[1] = 2.0 * r_tmp;
      } else {			/* beta_i *= (1+i), r_i *= (1+i) --> prod = 2i */
	beta_i[0] = beta_i_r;
	beta_i[1] = beta_i_r;
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i*=2i --> prod = 2i */
      r_i[0] = 0.0;
      r_i[1] = 2.0 * r_tmp;
    } else {			/* beta_i *= (1+i), r_i *= (1+i) */
      beta_i[0] = beta_i_r;
      beta_i[1] = beta_i_r;
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = beta_i_r;
    r_i[0] = r_tmp;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }
  for (i = 0; i < n; ++i) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[2 * i] = y_vec[i];
    y_i[2 * i + 1] = y_vec[i];
  }


  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_cdot2_s_c_testgen */

void BLAS_cdot2_c_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, float *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot2_c_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *r_i = (float *) r;
  float *head_x_i = (float *) head_x;
  float *tail_x_i = (float *) tail_x;
  float *y_i = y;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[2 * i];
    tail_x_vec[i] = tail_x_i[2 * i];
    y_vec[i] = y_i[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[2 * i];
    tail_x_vec[i] = tail_x_i[2 * i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i is 0 or 1 --> r_i *= 2i --> prod = 2i */
	r_i[0] = 0.0;
	r_i[1] = 2.0 * r_tmp;
      } else {			/* beta_i *= (1+i), r_i *= (1+i) --> prod = 2i */
	beta_i[0] = beta_i_r;
	beta_i[1] = beta_i_r;
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i*=2i --> prod = 2i */
      r_i[0] = 0.0;
      r_i[1] = 2.0 * r_tmp;
    } else {			/* beta_i *= (1+i), r_i *= (1+i) */
      beta_i[0] = beta_i_r;
      beta_i[1] = beta_i_r;
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = beta_i_r;
    r_i[0] = r_tmp;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }
  for (i = 0; i < n; ++i) {
    head_x_i[2 * i] = head_x_vec[i];
    head_x_i[2 * i + 1] = head_x_vec[i];
    tail_x_i[2 * i] = tail_x_vec[i];
    tail_x_i[2 * i + 1] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }
  if (conj == blas_conj) {
    for (i = 0; i < n; ++i) {
      head_x_i[2 * i + 1] = -head_x_i[2 * i + 1];
      tail_x_i[2 * i + 1] = -tail_x_i[2 * i + 1];
    }
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_cdot2_c_s_testgen */

void BLAS_zdot2_d_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    double *head_x, double *tail_x, double *y,
			    int *seed, void *r, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_d_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *head_x_i = head_x;
  double *tail_x_i = tail_x;
  double *y_i = y;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_tmp;
  double *head_x_vec;
  double *tail_x_vec;
  double *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (double *) blas_malloc(3 * n * sizeof(double));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 && ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.))) {	/* beta_i == 0 or 1 */
	r_i[0] = r_tmp;
	r_i[1] = 0.0;
      } else {			/* beta_i *= (1-i), r_i *= (1+i)/2 --> prod = 1 */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = r_tmp / 2.;
	r_i[1] = r_tmp / 2.;
      }
      r_true_l[1] = r_true_t[1] = 0.0;
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i *= (1+i) */
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    } else {			/* beta_i *= (1-i), r_i *= i */
      beta_i[0] = beta_i_r;
      beta_i[1] = -beta_i_r;
      r_i[0] = 0.;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = -beta_i_r;
    r_i[0] = 0;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  }
  for (i = 0; i < n; ++i) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }


  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_d_d_testgen */

void BLAS_zdot2_d_z_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    double *head_x, double *tail_x, void *y,
			    int *seed, void *r, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_d_z{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *head_x_i = head_x;
  double *tail_x_i = tail_x;
  double *y_i = (double *) y;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_tmp;
  double *head_x_vec;
  double *tail_x_vec;
  double *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (double *) blas_malloc(3 * n * sizeof(double));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[2 * i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i is 0 or 1 --> r_i *= 2i --> prod = 2i */
	r_i[0] = 0.0;
	r_i[1] = 2.0 * r_tmp;
      } else {			/* beta_i *= (1+i), r_i *= (1+i) --> prod = 2i */
	beta_i[0] = beta_i_r;
	beta_i[1] = beta_i_r;
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i*=2i --> prod = 2i */
      r_i[0] = 0.0;
      r_i[1] = 2.0 * r_tmp;
    } else {			/* beta_i *= (1+i), r_i *= (1+i) */
      beta_i[0] = beta_i_r;
      beta_i[1] = beta_i_r;
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = beta_i_r;
    r_i[0] = r_tmp;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }
  for (i = 0; i < n; ++i) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[2 * i] = y_vec[i];
    y_i[2 * i + 1] = y_vec[i];
  }


  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_d_z_testgen */

void BLAS_zdot2_z_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, double *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_z_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *head_x_i = (double *) head_x;
  double *tail_x_i = (double *) tail_x;
  double *y_i = y;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_tmp;
  double *head_x_vec;
  double *tail_x_vec;
  double *y_vec;

  alpha_i_r = alpha_i[0];
  alpha_i_i = alpha_i[1];
  beta_i_r = beta_i[0];
  beta_i_i = beta_i[1];
  head_x_vec = (double *) blas_malloc(3 * n * sizeof(double));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + n;
  y_vec = tail_x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    head_x_vec[i] = head_x_i[2 * i];
    tail_x_vec[i] = tail_x_i[2 * i];
    y_vec[i] = y_i[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    head_x_vec[i] = head_x_i[2 * i];
    tail_x_vec[i] = tail_x_i[2 * i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_i_r, alpha_flag,
		     &beta_i_r, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha_i is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha_i == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha_i == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i == 0 or 1 --> r_i *= (1+i) */
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      } else {			/* beta_i *= (1-i), r_i *= i --> prod = 1+i */
	beta_i[0] = beta_i_r;
	beta_i[1] = -beta_i_r;
	r_i[0] = 0.0;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha_i is a fixed multiple of (1+i) */
      alpha_i[0] = alpha_i_r;
      alpha_i[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta_i is 0 or 1 --> r_i *= 2i --> prod = 2i */
	r_i[0] = 0.0;
	r_i[1] = 2.0 * r_tmp;
      } else {			/* beta_i *= (1+i), r_i *= (1+i) --> prod = 2i */
	beta_i[0] = beta_i_r;
	beta_i[1] = beta_i_r;
	r_i[0] = r_tmp;
	r_i[1] = r_tmp;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r_i*=2i --> prod = 2i */
      r_i[0] = 0.0;
      r_i[1] = 2.0 * r_tmp;
    } else {			/* beta_i *= (1+i), r_i *= (1+i) */
      beta_i[0] = beta_i_r;
      beta_i[1] = beta_i_r;
      r_i[0] = r_tmp;
      r_i[1] = r_tmp;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha_i and beta_i are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    alpha_i[0] = alpha_i_r;
    alpha_i[1] = alpha_i_r;
    beta_i[0] = beta_i_r;
    beta_i[1] = beta_i_r;
    r_i[0] = r_tmp;
    r_i[1] = r_tmp;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }
  for (i = 0; i < n; ++i) {
    head_x_i[2 * i] = head_x_vec[i];
    head_x_i[2 * i + 1] = head_x_vec[i];
    tail_x_i[2 * i] = tail_x_vec[i];
    tail_x_i[2 * i + 1] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }
  if (conj == blas_conj) {
    for (i = 0; i < n; ++i) {
      head_x_i[2 * i + 1] = -head_x_i[2 * i + 1];
      tail_x_i[2 * i + 1] = -tail_x_i[2 * i + 1];
    }
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_z_d_testgen */

void BLAS_ddot2_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    float *head_x, float *tail_x, float *y, int *seed,
			    double *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot2_s_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *r_i = r;
  float *head_x_i = head_x;
  float *tail_x_i = tail_x;
  float *y_i = y;
  float alpha_tmp;
  float beta_tmp;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp = *alpha_i;
  beta_tmp = *beta_i;

  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_tmp, alpha_flag,
		     &beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     r_true_l, r_true_t);

  *alpha_i = alpha_tmp;
  *beta_i = beta_tmp;
  *r_i = r_tmp;
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_ddot2_s_s_testgen */

void BLAS_ddot2_s_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    float *head_x, float *tail_x, double *y,
			    int *seed, double *r, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot2_s_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *r_i = r;
  float *head_x_i = head_x;
  float *tail_x_i = tail_x;
  double *y_i = y;
  float alpha_tmp;
  float beta_tmp;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp = *alpha_i;
  beta_tmp = *beta_i;

  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_tmp, alpha_flag,
		     &beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     r_true_l, r_true_t);

  *alpha_i = alpha_tmp;
  *beta_i = beta_tmp;
  *r_i = r_tmp;
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_ddot2_s_d_testgen */

void BLAS_ddot2_d_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    double *head_x, double *tail_x, float *y,
			    int *seed, double *r, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot2_d_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = alpha;
  double *beta_i = beta;
  double *r_i = r;
  double *head_x_i = head_x;
  double *tail_x_i = tail_x;
  float *y_i = y;
  float alpha_tmp;
  float beta_tmp;
  float r_tmp;
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp = *alpha_i;
  beta_tmp = *beta_i;

  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float));
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
    y_vec[i] = y_i[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    tail_x_vec[i] = tail_x_i[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot2(n, n_fix2, n_mix, norm, conj,
		     &alpha_tmp, alpha_flag,
		     &beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, &r_tmp,
		     r_true_l, r_true_t);

  *alpha_i = alpha_tmp;
  *beta_i = beta_tmp;
  *r_i = r_tmp;
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    tail_x_i[i] = tail_x_vec[i];
    y_i[i] = y_vec[i];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_ddot2_d_s_testgen */

void BLAS_zdot2_c_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_c_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  float *head_x_i = (float *) head_x;
  float *tail_x_i = (float *) tail_x;
  float *y_i = (float *) y;
  float alpha_tmp[2];
  float beta_tmp[2];
  float r_tmp[2];
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp[0] = alpha_i[0];
  alpha_tmp[1] = alpha_i[1];
  beta_tmp[0] = beta_i[0];
  beta_tmp[1] = beta_i[1];
  inc *= 2;
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float) * 2);
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
    y_vec[i] = y_i[i];
    y_vec[i + 1] = y_i[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot2(n, n_fix2, n_mix, norm, conj,
		     alpha_tmp, alpha_flag,
		     beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, r_tmp,
		     r_true_l, r_true_t);

  alpha_i[0] = alpha_tmp[0];
  alpha_i[1] = alpha_tmp[1];
  beta_i[0] = beta_tmp[0];
  beta_i[1] = beta_tmp[1];
  r_i[0] = r_tmp[0];
  r_i[1] = r_tmp[1];
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    head_x_i[i + 1] = head_x_vec[i + 1];
    tail_x_i[i] = tail_x_vec[i];
    tail_x_i[i + 1] = tail_x_vec[i + 1];
    y_i[i] = y_vec[i];
    y_i[i + 1] = y_vec[i + 1];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_c_c_testgen */

void BLAS_zdot2_c_z_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_c_z{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  float *head_x_i = (float *) head_x;
  float *tail_x_i = (float *) tail_x;
  double *y_i = (double *) y;
  float alpha_tmp[2];
  float beta_tmp[2];
  float r_tmp[2];
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp[0] = alpha_i[0];
  alpha_tmp[1] = alpha_i[1];
  beta_tmp[0] = beta_i[0];
  beta_tmp[1] = beta_i[1];
  inc *= 2;
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float) * 2);
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
    y_vec[i] = y_i[i];
    y_vec[i + 1] = y_i[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot2(n, n_fix2, n_mix, norm, conj,
		     alpha_tmp, alpha_flag,
		     beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, r_tmp,
		     r_true_l, r_true_t);

  alpha_i[0] = alpha_tmp[0];
  alpha_i[1] = alpha_tmp[1];
  beta_i[0] = beta_tmp[0];
  beta_i[1] = beta_tmp[1];
  r_i[0] = r_tmp[0];
  r_i[1] = r_tmp[1];
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    head_x_i[i + 1] = head_x_vec[i + 1];
    tail_x_i[i] = tail_x_vec[i];
    tail_x_i[i + 1] = tail_x_vec[i + 1];
    y_i[i] = y_vec[i];
    y_i[i + 1] = y_vec[i + 1];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_c_z_testgen */

void BLAS_zdot2_z_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot2_z_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT2.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *head_x_i = (double *) head_x;
  double *tail_x_i = (double *) tail_x;
  float *y_i = (float *) y;
  float alpha_tmp[2];
  float beta_tmp[2];
  float r_tmp[2];
  float *head_x_vec;
  float *tail_x_vec;
  float *y_vec;

  alpha_tmp[0] = alpha_i[0];
  alpha_tmp[1] = alpha_i[1];
  beta_tmp[0] = beta_i[0];
  beta_tmp[1] = beta_i[1];
  inc *= 2;
  head_x_vec = (float *) blas_malloc(3 * n * sizeof(float) * 2);
  if (3 * n > 0 && head_x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_vec = head_x_vec + inc * n;
  y_vec = tail_x_vec + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
    y_vec[i] = y_i[i];
    y_vec[i + 1] = y_i[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    head_x_vec[i] = head_x_i[i];
    head_x_vec[i + 1] = head_x_i[i + 1];
    tail_x_vec[i] = tail_x_i[i];
    tail_x_vec[i + 1] = tail_x_i[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot2(n, n_fix2, n_mix, norm, conj,
		     alpha_tmp, alpha_flag,
		     beta_tmp, beta_flag,
		     head_x_vec, tail_x_vec, y_vec, seed, r_tmp,
		     r_true_l, r_true_t);

  alpha_i[0] = alpha_tmp[0];
  alpha_i[1] = alpha_tmp[1];
  beta_i[0] = beta_tmp[0];
  beta_i[1] = beta_tmp[1];
  r_i[0] = r_tmp[0];
  r_i[1] = r_tmp[1];
  for (i = 0; i < inc * n; i += inc) {
    head_x_i[i] = head_x_vec[i];
    head_x_i[i + 1] = head_x_vec[i + 1];
    tail_x_i[i] = tail_x_vec[i];
    tail_x_i[i + 1] = tail_x_vec[i + 1];
    y_i[i] = y_vec[i];
    y_i[i + 1] = y_vec[i + 1];
  }

  blas_free(head_x_vec);	/* also y_vec */
}				/* end BLAS_zdot2_z_c_testgen */
