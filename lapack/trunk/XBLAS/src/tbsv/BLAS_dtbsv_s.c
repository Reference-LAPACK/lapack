#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dtbsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, double alpha, const float *t, int ldt,
		  double *x, int incx)

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
 * alpha  (input) double
 * 
 * t      (input) float*
 *        Triangular Banded matrix
 *
 * x      (input) const double*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_dtbsv_s";

  int i, j;			/* used to keep track of loop counts */
  int xi;			/* used to index vector x */
  int start_xi;			/* used as the starting idx to vector x */
  int incxi;
  int Tij;			/* index inside of Banded structure */
  int dot_start, dot_start_inc1, dot_start_inc2, dot_inc;

  const float *t_i = t;		/* internal matrix t */
  double *x_i = x;		/* internal x */
  double alpha_i = alpha;	/* internal alpha */

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


  /* configuring the vector starting idx */
  if (incxi < 0) {
    start_xi = (1 - n) * incxi;
  } else {
    start_xi = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i == 0.0) {
    xi = start_xi;
    for (i = 0; i < n; i++) {
      x_i[xi] = 0.0;
      xi += incxi;
    }
    return;
  }
  /* check to see if k=0.  if so, we can optimize somewhat */
  if (k == 0) {
    if ((alpha_i == 1.0) && (diag == blas_unit_diag)) {
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







  {

    {
      double temp1;		/* temporary variable for calculations */
      double temp2;		/* temporary variable for calculations */
      double x_elem;
      float T_element;






      /*loop 1 */
      xi = start_xi;
      for (j = 0; j < k; j++) {

	/* each time through loop, xi lands on next x to compute. */
	x_elem = x_i[xi];
	/* preform the multiplication -
	   in this implementation we do not seperate the alpha = 1 case */
	temp1 = x_elem * alpha_i;

	xi = start_xi;

	Tij = dot_start;
	dot_start += dot_start_inc1;

	for (i = j; i > 0; i--) {
	  T_element = t_i[Tij];

	  x_elem = x_i[xi];
	  temp2 = x_elem * T_element;
	  temp1 = temp1 + (-temp2);
	  xi += incxi;
	  Tij += dot_inc;
	}			/* for across row */


	/* if the diagonal entry is not equal to one, then divide Xj by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = t_i[Tij];


	  temp1 = temp1 / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[xi] = temp1;
	xi += incxi;
      }				/* for j<k */
      /*end loop 1 */

      /*loop 2 continue without changing j to start */
      for (; j < n; j++) {

	/* each time through loop, xi lands on next x to compute. */
	x_elem = x_i[xi];
	temp1 = x_elem * alpha_i;

	xi = start_xi;
	start_xi += incxi;

	Tij = dot_start;
	dot_start += dot_start_inc2;

	for (i = k; i > 0; i--) {
	  T_element = t_i[Tij];

	  x_elem = x_i[xi];
	  temp2 = x_elem * T_element;
	  temp1 = temp1 + (-temp2);
	  xi += incxi;
	  Tij += dot_inc;
	}			/* for across row */


	/* if the diagonal entry is not equal to one, then divide by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = t_i[Tij];


	  temp1 = temp1 / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[xi] = temp1;
	xi += incxi;
      }				/* for j<n */


    }
  }
}				/* end BLAS_dtbsv_s */
