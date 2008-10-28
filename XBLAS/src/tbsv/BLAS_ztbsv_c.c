#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_ztbsv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const void *t, int ldt,
		  void *x, int incx)

/*
 * Purpose
 * =======
 * 
 * This routine solves :
 * 
 *     x <- alpha * inverse(t) * x
 * 
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major (blas_rowmajor, blas_colmajor)
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower (blas_upper, blas_lower)
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit (blas_unit_diag, blas_non_unit_diag)
 *
 * n      (input) int
 *        the dimension of t
 * 
 * k      (input) int
 *        the number of subdiagonals/superdiagonals of t
 *
 * alpha  (input) const void*
 * 
 * t      (input) void*
 *        Triangular Banded matrix
 *
 * x      (input) const void*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_ztbsv_c";

  int i, j;			/* used to keep track of loop counts */
  int xi;			/* used to index vector x */
  int start_xi;			/* used as the starting idx to vector x */
  int incxi;
  int Tij;			/* index inside of Banded structure */
  int dot_start, dot_start_inc1, dot_start_inc2, dot_inc;

  const float *t_i = (float *) t;	/* internal matrix t */
  double *x_i = (double *) x;	/* internal x */
  double *alpha_i = (double *) alpha;	/* internal alpha */

  if (order != blas_rowmajor && order != blas_colmajor) {
    BLAS_error(routine_name, -1, order, 0);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, 0);
  }
  if ((trans != blas_trans) && (trans != blas_no_trans) &&
      (trans != blas_conj) && (trans != blas_conj_trans)) {
    BLAS_error(routine_name, -2, uplo, 0);
  }
  if (diag != blas_non_unit_diag && diag != blas_unit_diag) {
    BLAS_error(routine_name, -4, diag, 0);
  }
  if (n < 0) {
    BLAS_error(routine_name, -5, n, 0);
  }
  if (k >= n) {
    BLAS_error(routine_name, -6, k, 0);
  }
  if ((ldt < 1) || (ldt <= k)) {
    BLAS_error(routine_name, -9, ldt, 0);
  }
  if (incx == 0) {
    BLAS_error(routine_name, -11, incx, 0);
  }

  if (n <= 0)
    return;

  incxi = incx;
  incxi *= 2;

  /* configuring the vector starting idx */
  if (incxi < 0) {
    start_xi = (1 - n) * incxi;
  } else {
    start_xi = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    xi = start_xi;
    for (i = 0; i < n; i++) {
      x_i[xi] = 0.0;
      x_i[xi + 1] = 0.0;
      xi += incxi;
    }
    return;
  }
  /* check to see if k=0.  if so, we can optimize somewhat */
  if (k == 0) {
    if (((alpha_i[0] == 1.0 && alpha_i[1] == 0.0))
	&& (diag == blas_unit_diag)) {
      /* nothing to do */
      return;
    } else {
      /* just run the loops as is. */

    }
  }

  /* get index variables prepared */
  if (((trans == blas_trans) || (trans == blas_conj_trans)) ^
      (order == blas_rowmajor)) {
    dot_start = k;
  } else {
    dot_start = 0;
  }

  if (((trans == blas_trans) || (trans == blas_conj_trans)) ^
      (order == blas_rowmajor)) {
    dot_inc = 1;
    dot_start_inc1 = ldt - 1;
    dot_start_inc2 = ldt;
  } else {
    dot_inc = ldt - 1;
    dot_start_inc1 = 1;
    dot_start_inc2 = ldt;
  }

  if (((trans == blas_trans) || (trans == blas_conj_trans)) ^
      (uplo == blas_lower)) {
    /*start at the first element of x */
    /* substitution will proceed forwards (forwardsubstitution) */
  } else {
    /*start at the last element of x */
    /* substitution will proceed backwards (backsubstitution) */
    dot_inc = -dot_inc;
    dot_start_inc1 = -dot_start_inc1;
    dot_start_inc2 = -dot_start_inc2;
    dot_start = ldt * (n - 1) + k - dot_start;
    /*order of the following 2 statements matters! */
    start_xi = start_xi + (n - 1) * incxi;
    incxi = -incxi;
  }

  dot_inc *= 2;
  dot_start *= 2;
  dot_start_inc1 *= 2;
  dot_start_inc2 *= 2;


  {

    {
      double temp1[2];		/* temporary variable for calculations */
      double temp2[2];		/* temporary variable for calculations */
      double x_elem[2];
      float T_element[2];




      if ((trans == blas_conj) || (trans == blas_conj_trans)) {
	/* conjugated */


	/*loop 1 */
	xi = start_xi;
	for (j = 0; j < k; j++) {

	  /* each time through loop, xi lands on next x to compute. */
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  /* preform the multiplication -
	     in this implementation we do not seperate the alpha = 1 case */
	  {
	    temp1[0] =
	      (double) x_elem[0] * alpha_i[0] -
	      (double) x_elem[1] * alpha_i[1];
	    temp1[1] =
	      (double) x_elem[0] * alpha_i[1] +
	      (double) x_elem[1] * alpha_i[0];
	  }

	  xi = start_xi;

	  Tij = dot_start;
	  dot_start += dot_start_inc1;

	  for (i = j; i > 0; i--) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];
	    T_element[1] = -T_element[1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      temp2[0] =
		(double) x_elem[0] * T_element[0] -
		(double) x_elem[1] * T_element[1];
	      temp2[1] =
		(double) x_elem[0] * T_element[1] +
		(double) x_elem[1] * T_element[0];
	    }
	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    xi += incxi;
	    Tij += dot_inc;
	  }			/* for across row */


	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un, eps1, ov1, un1;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      eps1 = pow(2.0, -53.0);
	      un1 = pow(2.0, -1022.0);
	      ov1 = 1.79769313486231571e+308;
	      /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
	      abs_a = fabs(temp1[0]);
	      abs_b = fabs(temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov1 / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		t = 2.0 / (eps1 * eps1);
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
	  x_i[xi] = temp1[0];
	  x_i[xi + 1] = temp1[1];
	  xi += incxi;
	}			/* for j<k */
	/*end loop 1 */

	/*loop 2 continue without changing j to start */
	for (; j < n; j++) {

	  /* each time through loop, xi lands on next x to compute. */
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  {
	    temp1[0] =
	      (double) x_elem[0] * alpha_i[0] -
	      (double) x_elem[1] * alpha_i[1];
	    temp1[1] =
	      (double) x_elem[0] * alpha_i[1] +
	      (double) x_elem[1] * alpha_i[0];
	  }

	  xi = start_xi;
	  start_xi += incxi;

	  Tij = dot_start;
	  dot_start += dot_start_inc2;

	  for (i = k; i > 0; i--) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];
	    T_element[1] = -T_element[1];
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      temp2[0] =
		(double) x_elem[0] * T_element[0] -
		(double) x_elem[1] * T_element[1];
	      temp2[1] =
		(double) x_elem[0] * T_element[1] +
		(double) x_elem[1] * T_element[0];
	    }
	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    xi += incxi;
	    Tij += dot_inc;
	  }			/* for across row */


	  /* if the diagonal entry is not equal to one, then divide by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un, eps1, ov1, un1;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      eps1 = pow(2.0, -53.0);
	      un1 = pow(2.0, -1022.0);
	      ov1 = 1.79769313486231571e+308;
	      /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
	      abs_a = fabs(temp1[0]);
	      abs_b = fabs(temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov1 / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		t = 2.0 / (eps1 * eps1);
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
	  x_i[xi] = temp1[0];
	  x_i[xi + 1] = temp1[1];
	  xi += incxi;
	}			/* for j<n */

      } else {
	/* not conjugated */


	/*loop 1 */
	xi = start_xi;
	for (j = 0; j < k; j++) {

	  /* each time through loop, xi lands on next x to compute. */
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  /* preform the multiplication -
	     in this implementation we do not seperate the alpha = 1 case */
	  {
	    temp1[0] =
	      (double) x_elem[0] * alpha_i[0] -
	      (double) x_elem[1] * alpha_i[1];
	    temp1[1] =
	      (double) x_elem[0] * alpha_i[1] +
	      (double) x_elem[1] * alpha_i[0];
	  }

	  xi = start_xi;

	  Tij = dot_start;
	  dot_start += dot_start_inc1;

	  for (i = j; i > 0; i--) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];

	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      temp2[0] =
		(double) x_elem[0] * T_element[0] -
		(double) x_elem[1] * T_element[1];
	      temp2[1] =
		(double) x_elem[0] * T_element[1] +
		(double) x_elem[1] * T_element[0];
	    }
	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    xi += incxi;
	    Tij += dot_inc;
	  }			/* for across row */


	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];


	    {
	      double S = 1.0, eps, ov, un, eps1, ov1, un1;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      eps1 = pow(2.0, -53.0);
	      un1 = pow(2.0, -1022.0);
	      ov1 = 1.79769313486231571e+308;
	      /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
	      abs_a = fabs(temp1[0]);
	      abs_b = fabs(temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov1 / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		t = 2.0 / (eps1 * eps1);
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
	  x_i[xi] = temp1[0];
	  x_i[xi + 1] = temp1[1];
	  xi += incxi;
	}			/* for j<k */
	/*end loop 1 */

	/*loop 2 continue without changing j to start */
	for (; j < n; j++) {

	  /* each time through loop, xi lands on next x to compute. */
	  x_elem[0] = x_i[xi];
	  x_elem[1] = x_i[xi + 1];
	  {
	    temp1[0] =
	      (double) x_elem[0] * alpha_i[0] -
	      (double) x_elem[1] * alpha_i[1];
	    temp1[1] =
	      (double) x_elem[0] * alpha_i[1] +
	      (double) x_elem[1] * alpha_i[0];
	  }

	  xi = start_xi;
	  start_xi += incxi;

	  Tij = dot_start;
	  dot_start += dot_start_inc2;

	  for (i = k; i > 0; i--) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];

	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      temp2[0] =
		(double) x_elem[0] * T_element[0] -
		(double) x_elem[1] * T_element[1];
	      temp2[1] =
		(double) x_elem[0] * T_element[1] +
		(double) x_elem[1] * T_element[0];
	    }
	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    xi += incxi;
	    Tij += dot_inc;
	  }			/* for across row */


	  /* if the diagonal entry is not equal to one, then divide by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = t_i[Tij];
	    T_element[1] = t_i[Tij + 1];


	    {
	      double S = 1.0, eps, ov, un, eps1, ov1, un1;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      eps1 = pow(2.0, -53.0);
	      un1 = pow(2.0, -1022.0);
	      ov1 = 1.79769313486231571e+308;
	      /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
	      abs_a = fabs(temp1[0]);
	      abs_b = fabs(temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov1 / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un1 / eps1 * 2) {	/* scale up a, b */
		t = 2.0 / (eps1 * eps1);
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
	  x_i[xi] = temp1[0];
	  x_i[xi + 1] = temp1[1];
	  xi += incxi;
	}			/* for j<n */

      }

    }
  }
}				/* end BLAS_ztbsv_c */
