#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ztrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *T, int ldt,
		  void *x, int incx, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine solve :
 * 
 *     x <- alpha * inverse(T) * x
 * 
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit
 *
 * n      (input) int
 *        the dimension of T
 * 
 * alpha  (input) const void*
 * 
 * T      (input) void*
 *        Triangular matrix
 *
 * x      (input) const void*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  char *routine_name = "BLAS_ztrsv";

  int i, j;			/* used to idx matrix */
  int ix, jx;			/* used to idx vector x */
  int start_x;			/* used as the starting idx to vector x */
  const double *T_i = (double *) T;	/* internal matrix T */
  double *x_i = (double *) x;	/* internal x */
  double *alpha_i = (double *) alpha;	/* internal alpha */
  double T_element[2];		/* temporary variable for an element of matrix A */
  int incT = 1;			/* internal ldt */

  if ((order != blas_rowmajor && order != blas_colmajor) ||
      (uplo != blas_upper && uplo != blas_lower) ||
      (trans != blas_trans && trans !=
       blas_no_trans && trans != blas_conj_trans) ||
      (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
      (ldt < n) || (incx == 0)) {
    BLAS_error(routine_name, 0, 0, NULL);
  }

  if (n <= 0)
    return;

  incT *= 2;
  incx *= 2;
  /* configuring the vector starting idx */
  if (incx <= 0) {
    start_x = -(n - 1) * incx;
  } else {
    start_x = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    ix = start_x;
    for (i = 0; i < n; i++) {
      x_i[ix] = 0.0;
      x_i[ix + 1] = 0.0;
      ix += incx;
    }
    return;
  }
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:

    {
      double temp1[2];		/* temporary variable for calculations */
      double temp2[2];		/* temporary variable for calculations */
      double temp3[2];		/* temporary variable for calculations */

      if ((order == blas_rowmajor &&
	   trans == blas_no_trans && uplo == blas_upper) ||
	  (order == blas_colmajor &&
	   trans != blas_no_trans && uplo == blas_lower)) {
	if (trans == blas_conj_trans) {

	  jx = start_x + (n - 1) * incx;
	  for (j = n - 1; j >= 0; j--) {

	    /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x + (n - 1) * incx;
	    for (i = n - 1; i >= j + 1; i--) {
	      T_element[0] = T_i[i * incT + j * ldt * incT];
	      T_element[1] = T_i[i * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];
	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix -= incx;
	    }			/* for j<n */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];

	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];

	    jx -= incx;
	  }			/* for j>=0 */
	} else {

	  jx = start_x + (n - 1) * incx;
	  for (j = n - 1; j >= 0; j--) {

	    /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x + (n - 1) * incx;
	    for (i = n - 1; i >= j + 1; i--) {
	      T_element[0] = T_i[i * incT + j * ldt * incT];
	      T_element[1] = T_i[i * incT + j * ldt * incT + 1];

	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix -= incx;
	    }			/* for j<n */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];

	    jx -= incx;
	  }			/* for j>=0 */
	}
      } else if ((order == blas_rowmajor &&
		  trans == blas_no_trans && uplo == blas_lower) ||
		 (order == blas_colmajor &&
		  trans != blas_no_trans && uplo == blas_upper)) {
	if (trans == blas_conj_trans) {

	  jx = start_x;
	  for (j = 0; j < n; j++) {

	    /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    /* multiply by alpha */
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x;
	    for (i = 0; i < j; i++) {
	      T_element[0] = T_i[i * incT + j * ldt * incT];
	      T_element[1] = T_i[i * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];
	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix += incx;
	    }			/* for i<j */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];

	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];
	    jx += incx;
	  }			/* for j<n */
	} else {

	  jx = start_x;
	  for (j = 0; j < n; j++) {

	    /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    /* multiply by alpha */
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x;
	    for (i = 0; i < j; i++) {
	      T_element[0] = T_i[i * incT + j * ldt * incT];
	      T_element[1] = T_i[i * incT + j * ldt * incT + 1];

	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix += incx;
	    }			/* for i<j */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];
	    jx += incx;
	  }			/* for j<n */
	}
      } else if ((order == blas_rowmajor &&
		  trans != blas_no_trans && uplo == blas_lower) ||
		 (order == blas_colmajor &&
		  trans == blas_no_trans && uplo == blas_upper)) {
	if (trans == blas_conj_trans) {

	  jx = start_x + (n - 1) * incx;
	  for (j = n - 1; j >= 0; j--) {

	    /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x + (n - 1) * incx;
	    for (i = n - 1; i >= j + 1; i--) {
	      T_element[0] = T_i[j * incT + i * ldt * incT];
	      T_element[1] = T_i[j * incT + i * ldt * incT + 1];
	      T_element[1] = -T_element[1];
	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix -= incx;
	    }			/* for j<n */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];

	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];

	    jx -= incx;
	  }			/* for j>=0 */
	} else {

	  jx = start_x + (n - 1) * incx;
	  for (j = n - 1; j >= 0; j--) {

	    /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x + (n - 1) * incx;
	    for (i = n - 1; i >= j + 1; i--) {
	      T_element[0] = T_i[j * incT + i * ldt * incT];
	      T_element[1] = T_i[j * incT + i * ldt * incT + 1];

	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix -= incx;
	    }			/* for j<n */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];

	    jx -= incx;
	  }			/* for j>=0 */
	}
      } else if ((order == blas_rowmajor &&
		  trans != blas_no_trans && uplo == blas_upper) ||
		 (order == blas_colmajor &&
		  trans == blas_no_trans && uplo == blas_lower)) {
	if (trans == blas_conj_trans) {

	  jx = start_x;
	  for (j = 0; j < n; j++) {

	    /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    /* multiply by alpha */
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x;
	    for (i = 0; i < j; i++) {
	      T_element[0] = T_i[j * incT + i * ldt * incT];
	      T_element[1] = T_i[j * incT + i * ldt * incT + 1];
	      T_element[1] = -T_element[1];
	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix += incx;
	    }			/* for i<j */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	      T_element[1] = -T_element[1];

	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];
	    jx += incx;
	  }			/* for j<n */
	} else {

	  jx = start_x;
	  for (j = 0; j < n; j++) {

	    /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	       i=j+1 to n-1           */
	    temp3[0] = x_i[jx];
	    temp3[1] = x_i[1 + jx];
	    /* multiply by alpha */
	    {
	      temp1[0] =
		(double) temp3[0] * alpha_i[0] -
		(double) temp3[1] * alpha_i[1];
	      temp1[1] =
		(double) temp3[0] * alpha_i[1] +
		(double) temp3[1] * alpha_i[0];
	    }

	    ix = start_x;
	    for (i = 0; i < j; i++) {
	      T_element[0] = T_i[j * incT + i * ldt * incT];
	      T_element[1] = T_i[j * incT + i * ldt * incT + 1];

	      temp3[0] = x_i[ix];
	      temp3[1] = x_i[1 + ix];
	      {
		temp2[0] =
		  (double) temp3[0] * T_element[0] -
		  (double) temp3[1] * T_element[1];
		temp2[1] =
		  (double) temp3[0] * T_element[1] +
		  (double) temp3[1] * T_element[0];
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      ix += incx;
	    }			/* for i<j */

	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element[0] = T_i[j * incT + j * ldt * incT];
	      T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	      {
		double S = 1.0, eps, ov, un;
		double abs_a, abs_b, abs_c, abs_d, ab, cd;
		double r;
		double t;
		double q[2];

		eps = pow(2.0, -24.0);
		un = pow(2.0, -126.0);
		ov = pow(2.0, 128.0) * (1 - eps);
		abs_a = fabs(temp1[0]);
		abs_b = fabs(temp1[1]);
		abs_c = fabs((double) T_element[0]);
		abs_d = fabs((double) T_element[1]);
		ab = MAX(abs_a, abs_b);
		cd = MAX(abs_c, abs_d);

		/* Scaling */
		if (ab > ov / 16) {	/* scale down a, b */
		  temp1[0] /= 16;
		  temp1[1] /= 16;
		  S = S * 16;
		}
		if (cd > ov / 16) {	/* scale down c, d */
		  T_element[0] /= 16;
		  T_element[1] /= 16;
		  S = S / 16;
		}
		if (ab < un / eps * 2) {	/* scale up a, b */
		  t = 2.0 / (eps * eps);
		  temp1[0] *= t;
		  temp1[1] *= t;
		  S = S / t;
		}
		if (cd < un / eps * 2) {	/* scale up c, d */
		  t = 2.0 / (eps * eps);
		  T_element[0] *= t;
		  T_element[1] *= t;
		  S = S * t;
		}

		/* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
		if (abs_c > abs_d) {
		  r = T_element[1] / T_element[0];
		  t = 1 / (T_element[0] + T_element[1] * r);
		  q[0] = (temp1[0] + temp1[1] * r) * t;
		  q[1] = (temp1[1] - temp1[0] * r) * t;
		} else {
		  r = T_element[0] / T_element[1];
		  t = 1 / (T_element[1] + T_element[0] * r);
		  q[0] = (temp1[1] + temp1[0] * r) * t;
		  q[1] = (-temp1[0] + temp1[1] * r) * t;
		}
		/* Scale back */
		temp1[0] = q[0] * S;
		temp1[1] = q[1] * S;
	      }

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[jx] = temp1[0];
	    x_i[jx + 1] = temp1[1];
	    jx += incx;
	  }			/* for j<n */
	}
      }
    }
    break;
  case blas_prec_extra:
    {
      FPU_FIX_DECL;
      FPU_FIX_START;
      {
	{
	  int inc_intx;		/* inc for intx */
	  double head_temp1[2], tail_temp1[2];	/* temporary variable for calculations */
	  double head_temp2[2], tail_temp2[2];	/* temporary variable for calculations */
	  double head_temp3[2], tail_temp3[2];	/* temporary variable for calculations */
	  double *head_intx, *tail_intx;
	  /* copy of x used for calculations */

	  /* allocate space for intx */
	  head_intx = (double *) blas_malloc(n * sizeof(double) * 2);
	  tail_intx = (double *) blas_malloc(n * sizeof(double) * 2);
	  if (n > 0 && (head_intx == NULL || tail_intx == NULL)) {
	    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
	  }

	  /* since intx is for internal usage, set it to 1 and then adjust
	     it if necessary */
	  inc_intx = 1;
	  inc_intx *= 2;

	  /* copy x to intx */
	  ix = start_x;
	  jx = 0;
	  for (i = 0; i < n; i++) {
	    head_temp1[0] = x_i[ix];
	    tail_temp1[0] = 0.0;
	    head_temp1[1] = x_i[1 + ix];
	    tail_temp1[1] = 0.0;
	    head_intx[jx] = head_temp1[0];
	    tail_intx[jx] = tail_temp1[0];
	    head_intx[1 + jx] = head_temp1[1];
	    tail_intx[1 + jx] = tail_temp1[1];
	    ix += incx;
	    jx += inc_intx;
	  }

	  if ((order == blas_rowmajor &&
	       trans == blas_no_trans && uplo == blas_upper) ||
	      (order == blas_colmajor &&
	       trans != blas_no_trans && uplo == blas_lower)) {
	    if (trans == blas_conj_trans) {

	      jx = (n - 1) * inc_intx;
	      for (j = n - 1; j >= 0; j--) {

		/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = (n - 1) * inc_intx;
		for (i = n - 1; i >= j + 1; i--) {
		  T_element[0] = T_i[i * incT + j * ldt * incT];
		  T_element[1] = T_i[i * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];
		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix -= inc_intx;
		}		/* for j<n */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];

		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];

		jx -= inc_intx;
	      }			/* for j>=0 */
	    } else {

	      jx = (n - 1) * inc_intx;
	      for (j = n - 1; j >= 0; j--) {

		/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = (n - 1) * inc_intx;
		for (i = n - 1; i >= j + 1; i--) {
		  T_element[0] = T_i[i * incT + j * ldt * incT];
		  T_element[1] = T_i[i * incT + j * ldt * incT + 1];

		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix -= inc_intx;
		}		/* for j<n */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];


		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];

		jx -= inc_intx;
	      }			/* for j>=0 */
	    }
	  } else if ((order == blas_rowmajor &&
		      trans == blas_no_trans && uplo == blas_lower) ||
		     (order == blas_colmajor &&
		      trans != blas_no_trans && uplo == blas_upper)) {
	    if (trans == blas_conj_trans) {

	      jx = 0;
	      for (j = 0; j < n; j++) {

		/* compute Xj = Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = 0;
		for (i = 0; i < j; i++) {
		  T_element[0] = T_i[i * incT + j * ldt * incT];
		  T_element[1] = T_i[i * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];
		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix += inc_intx;
		}		/* for i<j */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];

		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];
		jx += inc_intx;
	      }			/* for j<n */
	    } else {

	      jx = 0;
	      for (j = 0; j < n; j++) {

		/* compute Xj = Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = 0;
		for (i = 0; i < j; i++) {
		  T_element[0] = T_i[i * incT + j * ldt * incT];
		  T_element[1] = T_i[i * incT + j * ldt * incT + 1];

		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix += inc_intx;
		}		/* for i<j */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];


		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];
		jx += inc_intx;
	      }			/* for j<n */
	    }
	  } else if ((order == blas_rowmajor &&
		      trans != blas_no_trans && uplo == blas_lower) ||
		     (order == blas_colmajor &&
		      trans == blas_no_trans && uplo == blas_upper)) {
	    if (trans == blas_conj_trans) {

	      jx = (n - 1) * inc_intx;
	      for (j = n - 1; j >= 0; j--) {

		/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = (n - 1) * inc_intx;
		for (i = n - 1; i >= j + 1; i--) {
		  T_element[0] = T_i[j * incT + i * ldt * incT];
		  T_element[1] = T_i[j * incT + i * ldt * incT + 1];
		  T_element[1] = -T_element[1];
		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix -= inc_intx;
		}		/* for j<n */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];

		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];

		jx -= inc_intx;
	      }			/* for j>=0 */
	    } else {

	      jx = (n - 1) * inc_intx;
	      for (j = n - 1; j >= 0; j--) {

		/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = (n - 1) * inc_intx;
		for (i = n - 1; i >= j + 1; i--) {
		  T_element[0] = T_i[j * incT + i * ldt * incT];
		  T_element[1] = T_i[j * incT + i * ldt * incT + 1];

		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix -= inc_intx;
		}		/* for j<n */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];


		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];

		jx -= inc_intx;
	      }			/* for j>=0 */
	    }
	  } else if ((order == blas_rowmajor &&
		      trans != blas_no_trans && uplo == blas_upper) ||
		     (order == blas_colmajor &&
		      trans == blas_no_trans && uplo == blas_lower)) {
	    if (trans == blas_conj_trans) {

	      jx = 0;
	      for (j = 0; j < n; j++) {

		/* compute Xj = Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = 0;
		for (i = 0; i < j; i++) {
		  T_element[0] = T_i[j * incT + i * ldt * incT];
		  T_element[1] = T_i[j * incT + i * ldt * incT + 1];
		  T_element[1] = -T_element[1];
		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix += inc_intx;
		}		/* for i<j */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];
		  T_element[1] = -T_element[1];

		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];
		jx += inc_intx;
	      }			/* for j<n */
	    } else {

	      jx = 0;
	      for (j = 0; j < n; j++) {

		/* compute Xj = Xj - SUM Aij(or Aji) * Xi
		   i=j+1 to n-1           */
		head_temp3[0] = head_intx[jx];
		head_temp3[1] = head_intx[1 + jx];
		tail_temp3[0] = tail_intx[jx];
		tail_temp3[1] = tail_intx[1 + jx];
		/* multiply by alpha */
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_temp3[0];
		  tail_a0 = tail_temp3[0];
		  head_a1 = head_temp3[1];
		  tail_a1 = tail_temp3[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a0 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a1 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  head_t2 = -head_t2;
		  tail_t2 = -tail_t2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[0] = head_t1;
		  tail_temp1[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * alpha_i[0];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = alpha_i[1] * split;
		    b1 = con - alpha_i[1];
		    b1 = con - b1;
		    b2 = alpha_i[1] - b1;

		    c11 = head_a0 * alpha_i[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * alpha_i[1];
		    t1 = c11 + c2;
		    t2 = (c2 - (t1 - c11)) + c21;

		    head_t2 = t1 + t2;
		    tail_t2 = t2 - (head_t2 - t1);
		  }
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_t1 + head_t2;
		    bv = s1 - head_t1;
		    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_t1 + tail_t2;
		    bv = t1 - tail_t1;
		    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_t1 = t1 + t2;
		    tail_t1 = t2 - (head_t1 - t1);
		  }
		  head_temp1[1] = head_t1;
		  tail_temp1[1] = tail_t1;
		}


		ix = 0;
		for (i = 0; i < j; i++) {
		  T_element[0] = T_i[j * incT + i * ldt * incT];
		  T_element[1] = T_i[j * incT + i * ldt * incT + 1];

		  head_temp3[0] = head_intx[ix];
		  head_temp3[1] = head_intx[1 + ix];
		  tail_temp3[0] = tail_intx[ix];
		  tail_temp3[1] = tail_intx[1 + ix];
		  {
		    /* Compute complex-extra = complex-extra * complex-double. */
		    double head_a0, tail_a0;
		    double head_a1, tail_a1;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    head_a0 = head_temp3[0];
		    tail_a0 = tail_temp3[0];
		    head_a1 = head_temp3[1];
		    tail_a1 = tail_temp3[1];
		    /* real part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a0 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a1 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    head_t2 = -head_t2;
		    tail_t2 = -tail_t2;
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[0] = head_t1;
		    tail_temp2[0] = tail_t1;
		    /* imaginary part */
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a1 * split;
		      a11 = con - head_a1;
		      a11 = con - a11;
		      a21 = head_a1 - a11;
		      con = T_element[0] * split;
		      b1 = con - T_element[0];
		      b1 = con - b1;
		      b2 = T_element[0] - b1;

		      c11 = head_a1 * T_element[0];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a1 * T_element[0];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    {
		      /* Compute double-double = double-double * double. */
		      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_a0 * split;
		      a11 = con - head_a0;
		      a11 = con - a11;
		      a21 = head_a0 - a11;
		      con = T_element[1] * split;
		      b1 = con - T_element[1];
		      b1 = con - b1;
		      b2 = T_element[1] - b1;

		      c11 = head_a0 * T_element[1];
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_a0 * T_element[1];
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_t2 = t1 + t2;
		      tail_t2 = t2 - (head_t2 - t1);
		    }
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_t1 + head_t2;
		      bv = s1 - head_t1;
		      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_t1 + tail_t2;
		      bv = t1 - tail_t1;
		      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_t1 = t1 + t2;
		      tail_t1 = t2 - (head_t1 - t1);
		    }
		    head_temp2[1] = head_t1;
		    tail_temp2[1] = tail_t1;
		  }

		  {
		    double head_at, tail_at;
		    double head_bt, tail_bt;
		    double head_ct, tail_ct;

		    /* Real part */
		    head_at = head_temp1[0];
		    tail_at = tail_temp1[0];
		    head_bt = -head_temp2[0];
		    tail_bt = -tail_temp2[0];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[0] = head_ct;
		    tail_temp1[0] = tail_ct;
		    /* Imaginary part */
		    head_at = head_temp1[1];
		    tail_at = tail_temp1[1];
		    head_bt = -head_temp2[1];
		    tail_bt = -tail_temp2[1];
		    {
		      /* Compute double-double = double-double + double-double. */
		      double bv;
		      double s1, s2, t1, t2;

		      /* Add two hi words. */
		      s1 = head_at + head_bt;
		      bv = s1 - head_at;
		      s2 = ((head_bt - bv) + (head_at - (s1 - bv)));

		      /* Add two lo words. */
		      t1 = tail_at + tail_bt;
		      bv = t1 - tail_at;
		      t2 = ((tail_bt - bv) + (tail_at - (t1 - bv)));

		      s2 += t1;

		      /* Renormalize (s1, s2)  to  (t1, s2) */
		      t1 = s1 + s2;
		      s2 = s2 - (t1 - s1);

		      t2 += s2;

		      /* Renormalize (t1, t2)  */
		      head_ct = t1 + t2;
		      tail_ct = t2 - (head_ct - t1);
		    }
		    head_temp1[1] = head_ct;
		    tail_temp1[1] = tail_ct;
		  }
		  ix += inc_intx;
		}		/* for i<j */

		/* if the diagonal entry is not equal to one, then divide Xj by 
		   the entry */
		if (diag == blas_non_unit_diag) {
		  T_element[0] = T_i[j * incT + j * ldt * incT];
		  T_element[1] = T_i[j * incT + j * ldt * incT + 1];


		  {
		    double S = 1.0, eps, ov, un, eps1, ov1, un1;
		    double abs_a, abs_b, abs_c, abs_d, ab, cd;
		    double s;
		    double r;
		    double head_t, tail_t;
		    double head_t1, tail_t1;
		    double head_t2, tail_t2;
		    double head_q[2], tail_q[2];

		    eps = pow(2.0, -53.0);	/* double precision */
		    un = pow(2.0, -1022.0);
		    ov = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
		    eps1 = pow(2.0, -104.0);	/* extra precision */
		    un1 = pow(2.0, -1022.0);
		    ov1 = 1.79769313486231571e+308;
		    /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
		    abs_a = fabs(head_temp1[0]);
		    abs_b = fabs(head_temp1[1]);
		    abs_c = fabs((double) T_element[0]);
		    abs_d = fabs((double) T_element[1]);
		    ab = MAX(abs_a, abs_b);
		    cd = MAX(abs_c, abs_d);

		    /* Scaling */
		    if (ab > ov1 / 16) {	/* scale down a, b */
		      {
			double head_a, tail_a;
			double head_b, tail_b;
			head_a = head_temp1[0];
			tail_a = tail_temp1[0];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[0] = head_b;
			tail_temp1[0] = tail_b;
			head_a = head_temp1[1];
			tail_a = tail_temp1[1];
			{
			  /* Compute double-double = double-double / double,
			     using a Newton iteration scheme. */
			  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

			  /* Compute a DP approximation to the quotient. */
			  t1 = head_a / 16.0;

			  /* Split t1 and b into two parts with at most 26 bits each,
			     using the Dekker-Veltkamp method. */
			  con = t1 * split;
			  t11 = con - (con - t1);
			  t21 = t1 - t11;
			  con = 16.0 * split;
			  b1 = con - (con - 16.0);
			  b2 = 16.0 - b1;

			  /* Compute t1 * b using Dekker method. */
			  t12 = t1 * 16.0;
			  t22 =
			    (((t11 * b1 - t12) + t11 * b2) + t21 * b1) +
			    t21 * b2;

			  /* Compute dda - (t12, t22) using Knuth trick. */
			  t11 = head_a - t12;
			  e = t11 - head_a;
			  t21 =
			    ((-t12 - e) + (head_a - (t11 - e))) + tail_a -
			    t22;

			  /* Compute high-order word of (t11, t21) and divide by b. */
			  t2 = (t11 + t21) / 16.0;

			  /* The result is t1 + t2, after normalization. */
			  head_b = t1 + t2;
			  tail_b = t2 - (head_b - t1);
			}
			head_temp1[1] = head_b;
			tail_temp1[1] = tail_b;
		      }
		      S = S * 16;
		    }
		    if (cd > ov / 16) {	/* scale down c, d */
		      T_element[0] /= 16;
		      T_element[1] /= 16;
		      S = S / 16;
		    }
		    if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		      s = 2.0 / (eps1 * eps1);
		      {
			/* Compute complex-extra = complex-extra * real. */
			double head_a0, tail_a0;
			double head_a1, tail_a1;
			double head_t, tail_t;
			head_a0 = head_temp1[0];
			tail_a0 = tail_temp1[0];
			head_a1 = head_temp1[1];
			tail_a1 = tail_temp1[1];
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a0 * split;
			  a11 = con - head_a0;
			  a11 = con - a11;
			  a21 = head_a0 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a0 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a0 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[0] = head_t;
			tail_temp1[0] = tail_t;
			{
			  /* Compute double-double = double-double * double. */
			  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			  con = head_a1 * split;
			  a11 = con - head_a1;
			  a11 = con - a11;
			  a21 = head_a1 - a11;
			  con = s * split;
			  b1 = con - s;
			  b1 = con - b1;
			  b2 = s - b1;

			  c11 = head_a1 * s;
			  c21 =
			    (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			    a21 * b2;

			  c2 = tail_a1 * s;
			  t1 = c11 + c2;
			  t2 = (c2 - (t1 - c11)) + c21;

			  head_t = t1 + t2;
			  tail_t = t2 - (head_t - t1);
			}
			head_temp1[1] = head_t;
			tail_temp1[1] = tail_t;
		      }

		      S = S / s;
		    }
		    if (cd < un / eps * 2) {	/* scale up c, d */
		      s = 2.0 / (eps * eps);
		      T_element[0] *= s;
		      T_element[1] *= s;
		      S = S * s;
		    }

		    /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
		    if (abs_c > abs_d) {
		      r = T_element[1] / T_element[0];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[1] * split;
			b1 = con - T_element[1];
			b1 = con - b1;
			b2 = T_element[1] - b1;

			head_t = r * T_element[1];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[0];
			e = t1 - head_t;
			t2 =
			  ((T_element[0] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t2;
			tail_bt = -tail_t2;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t1 + head_bt;
			  bv = s1 - head_t1;
			  s2 = ((head_bt - bv) + (head_t1 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t1 + tail_bt;
			  bv = t1 - tail_t1;
			  t2 = ((tail_bt - bv) + (tail_t1 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* b - a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    } else {
		      r = T_element[0] / T_element[1];
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = r * split;
			a1 = con - r;
			a1 = con - a1;
			a2 = r - a1;
			con = T_element[0] * split;
			b1 = con - T_element[0];
			b1 = con - b1;
			b2 = T_element[0] - b1;

			head_t = r * T_element[0];
			tail_t =
			  (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double-double = double-double + double. */
			double e, t1, t2;

			/* Knuth trick. */
			t1 = head_t + T_element[1];
			e = t1 - head_t;
			t2 =
			  ((T_element[1] - e) + (head_t - (t1 - e))) + tail_t;

			/* The result is t1 + t2, after normalization. */
			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double + double-double. */
			double bv;
			double s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_t2 + head_t1;
			bv = s1 - head_t2;
			s2 = ((head_t1 - bv) + (head_t2 - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_t2 + tail_t1;
			bv = t1 - tail_t2;
			t2 = ((tail_t1 - bv) + (tail_t2 - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }		/* b + a*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[0] = head_t2;
		      tail_q[0] = tail_t2;
		      head_t1 = head_temp1[1];
		      tail_t1 = tail_temp1[1];	/* b */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_t1 * split;
			a11 = con - head_t1;
			a11 = con - a11;
			a21 = head_t1 - a11;
			con = r * split;
			b1 = con - r;
			b1 = con - b1;
			b2 = r - b1;

			c11 = head_t1 * r;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_t1 * r;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t2 = t1 + t2;
			tail_t2 = t2 - (head_t2 - t1);
		      }
		      head_t1 = head_temp1[0];
		      tail_t1 = tail_temp1[0];	/* a */
		      {
			double head_bt, tail_bt;
			head_bt = -head_t1;
			tail_bt = -tail_t1;
			{
			  /* Compute double-double = double-double + double-double. */
			  double bv;
			  double s1, s2, t1, t2;

			  /* Add two hi words. */
			  s1 = head_t2 + head_bt;
			  bv = s1 - head_t2;
			  s2 = ((head_bt - bv) + (head_t2 - (s1 - bv)));

			  /* Add two lo words. */
			  t1 = tail_t2 + tail_bt;
			  bv = t1 - tail_t2;
			  t2 = ((tail_bt - bv) + (tail_t2 - (t1 - bv)));

			  s2 += t1;

			  /* Renormalize (s1, s2)  to  (t1, s2) */
			  t1 = s1 + s2;
			  s2 = s2 - (t1 - s1);

			  t2 += s2;

			  /* Renormalize (t1, t2)  */
			  head_t2 = t1 + t2;
			  tail_t2 = t2 - (head_t2 - t1);
			}
		      }		/* -a + b*r */
		      {
			double q1, q2, q3;
			double a1, a2, b1, b2;
			double p1, p2, c;
			double s1, s2, v;
			double t1, t2;
			double r1, r2;
			double cona, conb;

			q1 = head_t2 / head_t;	/*  approximate quotient */

			/*  Compute  q1 * b  */
			cona = q1 * split;
			conb = head_t * split;
			a1 = cona - (cona - q1);
			b1 = conb - (conb - head_t);
			a2 = q1 - a1;
			b2 = head_t - b1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q1 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q1 * tail_t;


			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  a - (p1, p2)    */
			s1 = head_t2 - p1;
			v = s1 - head_t2;
			s2 = (head_t2 - (s1 - v)) - (p1 + v);

			t1 = tail_t2 - p2;
			v = t1 - tail_t2;
			t2 = (tail_t2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			r1 = t1 + t2;
			r2 = t2 - (r1 - t1);


			/*  Compute the next quotient. */
			q2 = r1 / head_t;


			/*  Compute residual   r1 - q2 * b          */
			cona = q2 * split;
			a1 = cona - (cona - q2);
			a2 = q2 - a1;

			/*  (p1, p2) is the product of high order terms. */
			p1 = q2 * head_t;
			p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;

			/*  Compute the low-order term */
			c = q2 * tail_t;

			/*  Compute  (s1, s2) = (p1, p2) + c */
			s1 = p1 + c;
			v = s1 - p1;
			s2 = ((c - v) + (p1 - (s1 - v))) + p2;

			/*  Renormalize. */
			p1 = s1 + s2;
			p2 = s2 - (p1 - s1);


			/*  Compute  (r1, r2) - (p1, p2)    */
			s1 = r1 - p1;
			v = s1 - r1;
			s2 = (r1 - (s1 - v)) - (p1 + v);

			t1 = r2 - p2;
			v = t1 - r2;
			t2 = (r2 - (t1 - v)) - (p2 + v);

			s2 += t1;
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;
			s1 = t1 + t2;

			/*  Compute the last correction. */
			q3 = s1 / head_t;

			/* Renormalize q1, q2, q3. */
			s1 = q2 + q3;
			s2 = q3 - (s1 - q2);

			head_t2 = q1 + s1;
			t1 = s1 - (head_t2 - q1);

			tail_t2 = s2 + t1;

		      }
		      head_q[1] = head_t2;
		      tail_q[1] = tail_t2;
		    }
		    /* Scale back */
		    if (S == 1.0) {
		      head_temp1[0] = head_q[0];
		      tail_temp1[0] = tail_q[0];
		      head_temp1[1] = head_q[1];
		      tail_temp1[1] = tail_q[1];
		    } else {
		      /* Compute complex-extra = complex-extra * real. */
		      double head_a0, tail_a0;
		      double head_a1, tail_a1;
		      double head_t, tail_t;
		      head_a0 = head_q[0];
		      tail_a0 = tail_q[0];
		      head_a1 = head_q[1];
		      tail_a1 = tail_q[1];
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a0 * split;
			a11 = con - head_a0;
			a11 = con - a11;
			a21 = head_a0 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a0 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[0] = head_t;
		      tail_temp1[0] = tail_t;
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = S * split;
			b1 = con - S;
			b1 = con - b1;
			b2 = S - b1;

			c11 = head_a1 * S;
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * S;
			t1 = c11 + c2;
			t2 = (c2 - (t1 - c11)) + c21;

			head_t = t1 + t2;
			tail_t = t2 - (head_t - t1);
		      }
		      head_temp1[1] = head_t;
		      tail_temp1[1] = tail_t;
		    }

		  }

		}
		/* if (diag == blas_non_unit_diag) */
		head_intx[jx] = head_temp1[0];
		tail_intx[jx] = tail_temp1[0];
		head_intx[1 + jx] = head_temp1[1];
		tail_intx[1 + jx] = tail_temp1[1];
		jx += inc_intx;
	      }			/* for j<n */
	    }
	  }

	  /* copy the final results from intx to x */
	  ix = start_x;
	  jx = 0;
	  for (i = 0; i < n; i++) {
	    head_temp1[0] = head_intx[jx];
	    head_temp1[1] = head_intx[1 + jx];
	    tail_temp1[0] = tail_intx[jx];
	    tail_temp1[1] = tail_intx[1 + jx];
	    x_i[ix] = head_temp1[0];
	    x_i[ix + 1] = head_temp1[1];
	    ix += incx;
	    jx += inc_intx;
	  }

	  blas_free(head_intx);
	  blas_free(tail_intx);
	}
      }
      FPU_FIX_STOP;
    }
    break;
  }
}
