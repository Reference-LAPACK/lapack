#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void BLAS_ssum_testgen(int n, int norm, float *x, int *seed,
		       double *sum_true_l, double *sum_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ssum{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vector X.
 *
 * norm    (input) int
 *         = -1 : the vector is scaled with norms near underflow.
 *         = 0  : the vector has norms of order 1.
 *         = 1  : the vector is scaled with norms near overflow.
 *
 * x       (output) float* contains generated test values
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * sum_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * sum_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  {
    int i;
    int xi, incxi = 1;
    float alpha;
    float beta;
    float r;
    float x_elem;
    float *tmp;

    float *x_i = x;

    alpha = 1.0;
    beta = 0.0;

    tmp = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && tmp == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    x_elem = 1.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      tmp[xi] = x_elem;
    }

    /* Call generator now. */
    testgen_BLAS_sdot(n, 0, n, norm, blas_conj,
		      &alpha, 1, &beta, 1,
		      tmp, x_i, seed, &r, sum_true_l, sum_true_t);

    blas_free(tmp);
  }
}
void BLAS_dsum_testgen(int n, int norm, double *x, int *seed,
		       double *sum_true_l, double *sum_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_dsum{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vector X.
 *
 * norm    (input) int
 *         = -1 : the vector is scaled with norms near underflow.
 *         = 0  : the vector has norms of order 1.
 *         = 1  : the vector is scaled with norms near overflow.
 *
 * x       (output) double* contains generated test values
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * sum_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * sum_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  {
    int i;
    int xi, incxi = 1;
    double alpha;
    double beta;
    double r;
    double x_elem;
    double *tmp;

    double *x_i = x;

    alpha = 1.0;
    beta = 0.0;

    tmp = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && tmp == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    x_elem = 1.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      tmp[xi] = x_elem;
    }

    /* Call generator now. */
    testgen_BLAS_ddot(n, 0, n, norm, blas_conj,
		      &alpha, 1, &beta, 1,
		      tmp, x_i, seed, &r, sum_true_l, sum_true_t);

    blas_free(tmp);
  }
}
void BLAS_csum_testgen(int n, int norm, void *x, int *seed,
		       double *sum_true_l, double *sum_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_csum{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vector X.
 *
 * norm    (input) int
 *         = -1 : the vector is scaled with norms near underflow.
 *         = 0  : the vector has norms of order 1.
 *         = 1  : the vector is scaled with norms near overflow.
 *
 * x       (output) void* contains generated test values
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * sum_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * sum_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  {
    int i;
    int xi, incxi = 1;
    float alpha[2];
    float beta[2];
    float r[2];
    float x_elem[2];
    float *tmp;

    float *x_i = (float *) x;

    alpha[0] = 1.0;
    alpha[1] = 0.0;
    beta[0] = beta[1] = 0.0;

    tmp = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && tmp == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    incxi *= 2;

    x_elem[0] = 1.0;
    x_elem[1] = 0.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      tmp[xi] = x_elem[0];
      tmp[xi + 1] = x_elem[1];
    }

    /* Call generator now. */
    testgen_BLAS_cdot(n, 0, n, norm, blas_conj,
		      alpha, 1, beta, 1,
		      tmp, x_i, seed, &r, sum_true_l, sum_true_t);

    blas_free(tmp);
  }
}
void BLAS_zsum_testgen(int n, int norm, void *x, int *seed,
		       double *sum_true_l, double *sum_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zsum{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vector X.
 *
 * norm    (input) int
 *         = -1 : the vector is scaled with norms near underflow.
 *         = 0  : the vector has norms of order 1.
 *         = 1  : the vector is scaled with norms near overflow.
 *
 * x       (output) void* contains generated test values
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * sum_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * sum_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  {
    int i;
    int xi, incxi = 1;
    double alpha[2];
    double beta[2];
    double r[2];
    double x_elem[2];
    double *tmp;

    double *x_i = (double *) x;

    alpha[0] = 1.0;
    alpha[1] = 0.0;
    beta[0] = beta[1] = 0.0;

    tmp = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && tmp == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    incxi *= 2;

    x_elem[0] = 1.0;
    x_elem[1] = 0.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      tmp[xi] = x_elem[0];
      tmp[xi + 1] = x_elem[1];
    }

    /* Call generator now. */
    testgen_BLAS_zdot(n, 0, n, norm, blas_conj,
		      alpha, 1, beta, 1,
		      tmp, x_i, seed, &r, sum_true_l, sum_true_t);

    blas_free(tmp);
  }
}
