#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


void test_BLAS_sdot2(int n, enum blas_conj_type conj, float alpha, float beta,
		     float rin, float rout, double r_true_l, double r_true_t,
		     float *x, int incx, float *head_y, float *tail_y,
		     int incy, double eps_int, double un_int,
		     double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) float
 *
 * beta    (input) float
 *
 * rin     (input) float
 *
 * rout    (input) float
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double
 *         The leading part of the truth.
 *
 * r_true_t (input) double
 *         The trailing part of the truth.
 *
 * x       (input) float*
 *
 * head_y
 * tail_y  (input) float*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U;
  double un_d, un_accurate, un_out;

  /* Set the starting position */
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
    S1 += fabs(x[ix]);
    S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
    ix += incx;
    iy += incy;
  }
  S *= fabs(alpha);
  S += fabs(beta * rin);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_S);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
  tmp1 = fabs((rout - r_true_l) - r_true_t);

  /* underflow */
  U = 2 * fabs(alpha) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;

  *test_ratio = tmp1 / ((n + 2) * (eps_int + eps_accurate) * S
			+ eps_out * fabs(r_true_l) + U);
}				/* end of test_BLAS_sdot */
void test_BLAS_ddot2(int n, enum blas_conj_type conj, double alpha,
		     double beta, double rin, double rout, double r_true_l,
		     double r_true_t, double *x, int incx, double *head_y,
		     double *tail_y, int incy, double eps_int, double un_int,
		     double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) double
 *
 * beta    (input) double
 *
 * rin     (input) double
 *
 * rout    (input) double
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double
 *         The leading part of the truth.
 *
 * r_true_t (input) double
 *         The trailing part of the truth.
 *
 * x       (input) double*
 *
 * head_y
 * tail_y  (input) double*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U;
  double un_d, un_accurate, un_out;

  /* Set the starting position */
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
    S1 += fabs(x[ix]);
    S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
    ix += incx;
    iy += incy;
  }
  S *= fabs(alpha);
  S += fabs(beta * rin);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp1 = fabs((rout - r_true_l) - r_true_t);

  /* underflow */
  U = 2 * fabs(alpha) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;

  *test_ratio = tmp1 / ((n + 2) * (eps_int + eps_accurate) * S
			+ eps_out * fabs(r_true_l) + U);
}				/* end of test_BLAS_ddot */
void test_BLAS_cdot2(int n, enum blas_conj_type conj, const void *alpha,
		     const void *beta, const void *rin, const void *rout,
		     double *r_true_l, double *r_true_t, void *x, int incx,
		     void *head_y, void *tail_y, int incy, double eps_int,
		     double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = (float *) x;
  float *head_y_i = (float *) head_y;
  float *tail_y_i = (float *) tail_y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *rin_i = (float *) rin;
  float *rout_i = (float *) rout;
  float x_ii[2];
  float y_ii[2];

  /* Set the starting position */
  incx *= 2;
  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = x_ii[0] * y_ii[0] - x_ii[1] * y_ii[1];
      prod[1] = x_ii[0] * y_ii[1] + x_ii[1] * y_ii[0];
    }

    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = x_ii[0] * y_ii[0] - x_ii[1] * y_ii[1];
      prod[1] = x_ii[0] * y_ii[1] + x_ii[1] * y_ii[0];
    }
    /* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = beta_i[0] * rin_i[0] - beta_i[1] * rin_i[1];
    prod[1] = beta_i[0] * rin_i[1] + beta_i[1] * rin_i[0];
  }

  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_S);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_cdot */
void test_BLAS_zdot2(int n, enum blas_conj_type conj, const void *alpha,
		     const void *beta, const void *rin, const void *rout,
		     double *r_true_l, double *r_true_t, void *x, int incx,
		     void *head_y, void *tail_y, int incy, double eps_int,
		     double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  double *x_i = (double *) x;
  double *head_y_i = (double *) head_y;
  double *tail_y_i = (double *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  double x_ii[2];
  double y_ii[2];

  /* Set the starting position */
  incx *= 2;
  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_cdot2_s_s(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, float *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) float*
 *
 * head_y
 * tail_y  (input) float*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = x;
  float *head_y_i = head_y;
  float *tail_y_i = tail_y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *rin_i = (float *) rin;
  float *rout_i = (float *) rout;
  float x_ii;
  float y_ii;

  /* Set the starting position */


  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = head_y_i[iy];
    S1 += fabs(x_ii);
    S2 += fabs(y_ii);

    prod[0] = x_ii * y_ii;
    prod[1] = 0.0;
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii = tail_y_i[iy];
    prod[0] = x_ii * y_ii;
    prod[1] = 0.0;		/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += fabs(y_ii);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = beta_i[0] * rin_i[0] - beta_i[1] * rin_i[1];
    prod[1] = beta_i[0] * rin_i[1] + beta_i[1] * rin_i[0];
  }

  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_S);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_cdot */
void test_BLAS_cdot2_s_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, float *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) float*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = x;
  float *head_y_i = (float *) head_y;
  float *tail_y_i = (float *) tail_y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *rin_i = (float *) rin;
  float *rout_i = (float *) rout;
  float x_ii;
  float y_ii[2];

  /* Set the starting position */

  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    S1 += fabs(x_ii);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = y_ii[0] * x_ii;
      prod[1] = y_ii[1] * x_ii;
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = y_ii[0] * x_ii;
      prod[1] = y_ii[1] * x_ii;
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = beta_i[0] * rin_i[0] - beta_i[1] * rin_i[1];
    prod[1] = beta_i[0] * rin_i[1] + beta_i[1] * rin_i[0];
  }

  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_S);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_cdot */
void test_BLAS_cdot2_c_s(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) float*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = (float *) x;
  float *head_y_i = head_y;
  float *tail_y_i = tail_y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *rin_i = (float *) rin;
  float *rout_i = (float *) rout;
  float x_ii[2];
  float y_ii;

  /* Set the starting position */
  incx *= 2;

  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii = head_y_i[iy];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += fabs(y_ii);

    {
      prod[0] = x_ii[0] * y_ii;
      prod[1] = x_ii[1] * y_ii;
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii = tail_y_i[iy];
    {
      prod[0] = x_ii[0] * y_ii;
      prod[1] = x_ii[1] * y_ii;
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += fabs(y_ii);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = beta_i[0] * rin_i[0] - beta_i[1] * rin_i[1];
    prod[1] = beta_i[0] * rin_i[1] + beta_i[1] * rin_i[0];
  }

  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_S);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_cdot */
void test_BLAS_zdot2_d_d(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, double *x,
			 int incx, double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) double*
 *
 * head_y
 * tail_y  (input) double*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  double *x_i = x;
  double *head_y_i = head_y;
  double *tail_y_i = tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  double x_ii;
  double y_ii;

  /* Set the starting position */


  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = head_y_i[iy];
    S1 += fabs(x_ii);
    S2 += fabs(y_ii);

    prod[0] = x_ii * y_ii;
    prod[1] = 0.0;
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii = tail_y_i[iy];
    prod[0] = x_ii * y_ii;
    prod[1] = 0.0;		/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += fabs(y_ii);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_zdot2_d_z(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, double *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) double*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  double *x_i = x;
  double *head_y_i = (double *) head_y;
  double *tail_y_i = (double *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  double x_ii;
  double y_ii[2];

  /* Set the starting position */

  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    S1 += fabs(x_ii);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = y_ii[0] * x_ii;
      prod[1] = y_ii[1] * x_ii;
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = y_ii[0] * x_ii;
      prod[1] = y_ii[1] * x_ii;
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_zdot2_z_d(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) double*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  double *x_i = (double *) x;
  double *head_y_i = head_y;
  double *tail_y_i = tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  double x_ii[2];
  double y_ii;

  /* Set the starting position */
  incx *= 2;

  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii = head_y_i[iy];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += fabs(y_ii);

    {
      prod[0] = x_ii[0] * y_ii;
      prod[1] = x_ii[1] * y_ii;
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii = tail_y_i[iy];
    {
      prod[0] = x_ii[0] * y_ii;
      prod[1] = x_ii[1] * y_ii;
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += fabs(y_ii);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_ddot2_s_s(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, float *x, int incx,
			 float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) double
 *
 * beta    (input) double
 *
 * rin     (input) double
 *
 * rout    (input) double
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double
 *         The leading part of the truth.
 *
 * r_true_t (input) double
 *         The trailing part of the truth.
 *
 * x       (input) float*
 *
 * head_y
 * tail_y  (input) float*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U;
  double un_d, un_accurate, un_out;

  /* Set the starting position */
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
    S1 += fabs(x[ix]);
    S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
    ix += incx;
    iy += incy;
  }
  S *= fabs(alpha);
  S += fabs(beta * rin);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp1 = fabs((rout - r_true_l) - r_true_t);

  /* underflow */
  U = 2 * fabs(alpha) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;

  *test_ratio = tmp1 / ((n + 2) * (eps_int + eps_accurate) * S
			+ eps_out * fabs(r_true_l) + U);
}				/* end of test_BLAS_ddot */
void test_BLAS_ddot2_s_d(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, float *x, int incx,
			 double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) double
 *
 * beta    (input) double
 *
 * rin     (input) double
 *
 * rout    (input) double
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double
 *         The leading part of the truth.
 *
 * r_true_t (input) double
 *         The trailing part of the truth.
 *
 * x       (input) float*
 *
 * head_y
 * tail_y  (input) double*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U;
  double un_d, un_accurate, un_out;

  /* Set the starting position */
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
    S1 += fabs(x[ix]);
    S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
    ix += incx;
    iy += incy;
  }
  S *= fabs(alpha);
  S += fabs(beta * rin);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp1 = fabs((rout - r_true_l) - r_true_t);

  /* underflow */
  U = 2 * fabs(alpha) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;

  *test_ratio = tmp1 / ((n + 2) * (eps_int + eps_accurate) * S
			+ eps_out * fabs(r_true_l) + U);
}				/* end of test_BLAS_ddot */
void test_BLAS_ddot2_d_s(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, double *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) double
 *
 * beta    (input) double
 *
 * rin     (input) double
 *
 * rout    (input) double
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double
 *         The leading part of the truth.
 *
 * r_true_t (input) double
 *         The trailing part of the truth.
 *
 * x       (input) double*
 *
 * head_y
 * tail_y  (input) float*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U;
  double un_d, un_accurate, un_out;

  /* Set the starting position */
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
    S1 += fabs(x[ix]);
    S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
    ix += incx;
    iy += incy;
  }
  S *= fabs(alpha);
  S += fabs(beta * rin);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp1 = fabs((rout - r_true_l) - r_true_t);

  /* underflow */
  U = 2 * fabs(alpha) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;

  *test_ratio = tmp1 / ((n + 2) * (eps_int + eps_accurate) * S
			+ eps_out * fabs(r_true_l) + U);
}				/* end of test_BLAS_ddot */
void test_BLAS_zdot2_c_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = (float *) x;
  float *head_y_i = (float *) head_y;
  float *tail_y_i = (float *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  float x_ii[2];
  float y_ii[2];

  /* Set the starting position */
  incx *= 2;
  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_zdot2_c_z(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  float *x_i = (float *) x;
  double *head_y_i = (double *) head_y;
  double *tail_y_i = (double *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  float x_ii[2];
  double y_ii[2];

  /* Set the starting position */
  incx *= 2;
  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
void test_BLAS_zdot2_z_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio)

/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) const void*
 *
 * beta    (input) const void*
 *
 * rin     (input) const void*
 *
 * rout    (input) const void*
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) double*
 *         The leading part of the truth.
 *
 * r_true_t (input) double*
 *         The trailing part of the truth.
 *
 * x       (input) void*
 *
 * head_y
 * tail_y  (input) void*
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */
{
  int i, ix, iy;
  double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
  double un_d, un_accurate, un_out;
  double *x_i = (double *) x;
  float *head_y_i = (float *) head_y;
  float *tail_y_i = (float *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *rin_i = (double *) rin;
  double *rout_i = (double *) rout;
  double x_ii[2];
  float y_ii[2];

  /* Set the starting position */
  incx *= 2;
  incy *= 2;
  ix = 0;
  iy = 0;
  if (incx < 0)
    ix = -(n - 1) * incx;
  if (incy < 0)
    iy = -(n - 1) * incy;

  /* computing S */
  S = S1 = S2 = 0.;
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    y_ii[0] = head_y_i[iy];
    y_ii[1] = head_y_i[iy + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    S1 += sqrt(x_ii[0] * x_ii[0] + x_ii[1] * x_ii[1]);
    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);

    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    /* Now get the tail of y */
    y_ii[0] = tail_y_i[iy];
    y_ii[1] = tail_y_i[iy + 1];
    {
      prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
      prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
    }				/* prod = x[i]*y[i] */
    S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

    S2 += sqrt(y_ii[0] * y_ii[0] + y_ii[1] * y_ii[1]);
    ix += incx;
    iy += incy;
  }
  S *= sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]);
  {
    prod[0] = (double) beta_i[0] * rin_i[0] - (double) beta_i[1] * rin_i[1];
    prod[1] = (double) beta_i[0] * rin_i[1] + (double) beta_i[1] * rin_i[0];
  }
  S += sqrt(prod[0] * prod[0] + prod[1] * prod[1]);

  un_d = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	     (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  S = MAX(S, un_d);

  eps_accurate = power(2, -BITS_E);
  un_accurate = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		    (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
  eps_out = power(2, -BITS_D);
  un_out = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
	       (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
  tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
  tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
  tmp1 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

  /* underflow */
  U = 2 * sqrt(alpha_i[0] * alpha_i[0] + alpha_i[1] * alpha_i[1]) * n + 3;
  U = MAX(U, S1 + 2 * n + 1);
  U = MAX(U, S2 + 2 * n + 1) * (un_int + un_accurate) + un_out;
  U *= 2 * sqrt(2.);

  *test_ratio = tmp1 / (2 * sqrt(2.) * (n + 2) * (eps_int + eps_accurate) * S
			+
			sqrt(2.) * eps_out * sqrt(r_true_l[0] * r_true_l[0] +
						  r_true_l[1] * r_true_l[1]) +
			U);
}

  /* end of test_BLAS_zdot */
