#include <stdlib.h>
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void
BLAS_sdot_x_testgen(int n, int n_fix2, int n_mix, int norm,
                    enum blas_conj_type conj,
                    float *alpha, int alpha_flag,
                    float *beta, int beta_flag,
                    double *x_l, double *x_t,
                    float *y, int *seed, float *r,
                    double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates one row of the test inputs to BLAS_{s,d}_trsv_x.
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
 * x_l     (input) double*
 *         The leading part of the double-double X vector.
 *
 * x_t     (input) double*
 *         The trailing part of the double-double X vector.
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * r       (output) float*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  /* Call generator now. */
  testgen_BLAS_sdot_x(n, n_fix2, n_mix, norm, conj,
                      alpha, alpha_flag, beta, beta_flag,
                      x_l, x_t, y, seed, r, r_true_l, r_true_t);
}                               /* end BLAS_sdot_x_testgen */


void
BLAS_ddot_x_testgen(int n, int n_fix2, int n_mix, int norm,
                    enum blas_conj_type conj,
                    double *alpha, int alpha_flag,
                    double *beta, int beta_flag,
                    double *x_l, double *x_t,
                    double *y, int *seed, double *r,
                    double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test inputs to BLAS_ddot_s_d{_x}.
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
 * x_l     (input) double*
 *         The leading part of the X vector in double-double.
 *
 * x_t     (input) double*
 *         The trailing part of the X vector in double-double.
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
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
  float alpha_i;
  float beta_i;
  float r_i;
  float *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  y_i = (float *) malloc(n * sizeof(float));
  if (y_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < n; ++i)
    y_i[i] = y[i];

  /*printf("dot_trsv0> r=%e, r_i=%e\n", *r, r_i); */
  /* Call generator now. */
  testgen_BLAS_sdot_x(n, n_fix2, n_mix, norm, conj,
                      &alpha_i, alpha_flag, &beta_i, beta_flag,
                      x_l, x_t, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  /*printf("dot_trsv> r=%e, r_i=%e\n", *r, r_i); */

  for (i = 0; i < n; ++i)
    y[i] = y_i[i];

  free(y_i);
}                               /* end BLAS_ddot_x_testgen */


void
BLAS_ddot_s_x_testgen(int n, int n_fix2, int n_mix, int norm,
                      enum blas_conj_type conj,
                      double *alpha, int alpha_flag,
                      double *beta, int beta_flag,
                      double *x_l, double *x_t,
                      float *y, int *seed, double *r,
                      double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test inputs to BLAS_ddot_s_d{_x}.
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
 * x_l     (input) double*
 *         The leading part of the X vector in double-double.
 *
 * x_t     (input) double*
 *         The trailing part of the X vector in double-double.
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  float alpha_i;
  float beta_i;
  float r_i;

  alpha_i = *alpha;
  beta_i = *beta;

  /* Call generator now. */
  testgen_BLAS_sdot_x(n, n_fix2, n_mix, norm, conj,
                      &alpha_i, alpha_flag, &beta_i, beta_flag,
                      x_l, x_t, y, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
}                               /* end BLAS_ddot_s_x_testgen */
