#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void test_BLAS_ssum(int n, float sum_comp, double sum_true_l,
		    double sum_true_t, float *x, int incx,
		    double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from sum over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vector X.
 *
 * sum_comp    (input) float
 *             This result was computed by some other routine, and will be
 *             tested by this routine by comparing it with the truth.
 *
 * sum_true_l  (input) double
 *             leading part of true sum value
 *
 * sum_true_t  (input) double 
 *             tailing part of true sum value
 *
 * x       (input) float*
 *
 * incx    (input) int
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) double*
 *         The ratio of computed error for r over the error bound.
 */
{
  {
    int i;
    int xi;
    int incxi = 1;
    float y_elem;
    float alpha;
    float beta;
    float r;

    float *y;
    y = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && y == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    /* set each element of y to be 1.0 */
    y_elem = 1.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      y[xi] = y_elem;
    }

    alpha = 1.0;
    beta = 0.0;
    r = 0.0;

    test_BLAS_sdot(n, blas_no_conj, alpha, beta, r, sum_comp,
		   sum_true_l, sum_true_t, x, incx, y, 1,
		   eps_int, un_int, test_ratio);
    blas_free(y);
  }
}
void test_BLAS_dsum(int n, double sum_comp, double sum_true_l,
		    double sum_true_t, double *x, int incx,
		    double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from sum over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vector X.
 *
 * sum_comp    (input) double
 *             This result was computed by some other routine, and will be
 *             tested by this routine by comparing it with the truth.
 *
 * sum_true_l  (input) double
 *             leading part of true sum value
 *
 * sum_true_t  (input) double 
 *             tailing part of true sum value
 *
 * x       (input) double*
 *
 * incx    (input) int
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) double*
 *         The ratio of computed error for r over the error bound.
 */
{
  {
    int i;
    int xi;
    int incxi = 1;
    double y_elem;
    double alpha;
    double beta;
    double r;

    double *y;
    y = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && y == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }



    /* set each element of y to be 1.0 */
    y_elem = 1.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      y[xi] = y_elem;
    }

    alpha = 1.0;
    beta = 0.0;
    r = 0.0;

    test_BLAS_ddot(n, blas_no_conj, alpha, beta, r, sum_comp,
		   sum_true_l, sum_true_t, x, incx, y, 1,
		   eps_int, un_int, test_ratio);
    blas_free(y);
  }
}
void test_BLAS_csum(int n, const void *sum_comp, double *sum_true_l,
		    double *sum_true_t, void *x, int incx,
		    double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from sum over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vector X.
 *
 * sum_comp    (input) const void*
 *             This result was computed by some other routine, and will be
 *             tested by this routine by comparing it with the truth.
 *
 * sum_true_l  (input) double*
 *             leading part of true sum value
 *
 * sum_true_t  (input) double* 
 *             tailing part of true sum value
 *
 * x       (input) void*
 *
 * incx    (input) int
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) double*
 *         The ratio of computed error for r over the error bound.
 */
{
  {
    int i;
    int xi;
    int incxi = 1;
    float y_elem[2];
    float alpha[2];
    float beta[2];
    float r[2];

    float *y;
    y = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && y == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    incxi *= 2;

    /* set each element of y to be 1.0 */
    y_elem[0] = 1.0;
    y_elem[1] = 0.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      y[xi] = y_elem[0];
      y[xi + 1] = y_elem[1];
    }

    alpha[0] = 1.0;
    alpha[1] = 0.0;
    beta[0] = beta[1] = 0.0;
    r[0] = r[1] = 0.0;

    test_BLAS_cdot(n, blas_no_conj, alpha, beta, r, sum_comp,
		   sum_true_l, sum_true_t, x, incx, y, 1,
		   eps_int, un_int, test_ratio);
    blas_free(y);
  }
}
void test_BLAS_zsum(int n, const void *sum_comp, double *sum_true_l,
		    double *sum_true_t, void *x, int incx,
		    double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from sum over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vector X.
 *
 * sum_comp    (input) const void*
 *             This result was computed by some other routine, and will be
 *             tested by this routine by comparing it with the truth.
 *
 * sum_true_l  (input) double*
 *             leading part of true sum value
 *
 * sum_true_t  (input) double* 
 *             tailing part of true sum value
 *
 * x       (input) void*
 *
 * incx    (input) int
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) double*
 *         The ratio of computed error for r over the error bound.
 */
{
  {
    int i;
    int xi;
    int incxi = 1;
    double y_elem[2];
    double alpha[2];
    double beta[2];
    double r[2];

    double *y;
    y = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && y == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    incxi *= 2;

    /* set each element of y to be 1.0 */
    y_elem[0] = 1.0;
    y_elem[1] = 0.0;
    for (i = 0, xi = 0; i < n; i++, xi += incxi) {
      y[xi] = y_elem[0];
      y[xi + 1] = y_elem[1];
    }

    alpha[0] = 1.0;
    alpha[1] = 0.0;
    beta[0] = beta[1] = 0.0;
    r[0] = r[1] = 0.0;

    test_BLAS_zdot(n, blas_no_conj, alpha, beta, r, sum_comp,
		   sum_true_l, sum_true_t, x, incx, y, 1,
		   eps_int, un_int, test_ratio);
    blas_free(y);
  }
}
