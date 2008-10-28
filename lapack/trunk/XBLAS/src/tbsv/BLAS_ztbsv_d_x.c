#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_ztbsv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, int k, const void *alpha, const double *t, int ldt,
		    void *x, int incx, enum blas_prec_type prec)

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
 * t      (input) double*
 *        Triangular Banded matrix
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
  /* Routine name */
  static const char routine_name[] = "BLAS_ztbsv_d_x";

  int i, j;			/* used to keep track of loop counts */
  int xi;			/* used to index vector x */
  int start_xi;			/* used as the starting idx to vector x */
  int incxi;
  int Tij;			/* index inside of Banded structure */
  int dot_start, dot_start_inc1, dot_start_inc2, dot_inc;

  const double *t_i = t;	/* internal matrix t */
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
      /* must set prec to output. Ignore user input of prec */
      prec = blas_prec_double;
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







  switch (prec) {

  case blas_prec_single:
  case blas_prec_indigenous:
  case blas_prec_double:{
      {

	{
	  double temp1[2];	/* temporary variable for calculations */
	  double temp2[2];	/* temporary variable for calculations */
	  double x_elem[2];
	  double T_element;






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
	      T_element = t_i[Tij];

	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		temp2[0] = x_elem[0] * T_element;
		temp2[1] = x_elem[1] * T_element;
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      xi += incxi;
	      Tij += dot_inc;
	    }			/* for across row */


	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element = t_i[Tij];


	      temp1[0] = temp1[0] / T_element;
	      temp1[1] = temp1[1] / T_element;

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
	      T_element = t_i[Tij];

	      x_elem[0] = x_i[xi];
	      x_elem[1] = x_i[xi + 1];
	      {
		temp2[0] = x_elem[0] * T_element;
		temp2[1] = x_elem[1] * T_element;
	      }
	      temp1[0] = temp1[0] - temp2[0];
	      temp1[1] = temp1[1] - temp2[1];
	      xi += incxi;
	      Tij += dot_inc;
	    }			/* for across row */


	    /* if the diagonal entry is not equal to one, then divide by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element = t_i[Tij];


	      temp1[0] = temp1[0] / T_element;
	      temp1[1] = temp1[1] / T_element;

	    }
	    /* if (diag == blas_non_unit_diag) */
	    x_i[xi] = temp1[0];
	    x_i[xi + 1] = temp1[1];
	    xi += incxi;
	  }			/* for j<n */


	}
      }
      break;
    }

  case blas_prec_extra:{
      {

	{
	  double head_temp1[2], tail_temp1[2];	/* temporary variable for calculations */
	  double head_temp2[2], tail_temp2[2];	/* temporary variable for calculations */
	  double head_temp3[2], tail_temp3[2];	/* temporary variable for calculations */
	  double x_elem[2];
	  double T_element;	/* temporary variable for an element of matrix T */

	  int x_inti = 0, inc_x_inti = 1;
	  int k_compare = k;	/*used for comparisons with x_inti */
	  double *head_x_internal, *tail_x_internal;

	  FPU_FIX_DECL;

	  k_compare *= 2;
	  inc_x_inti *= 2;
	  head_x_internal = (double *) blas_malloc(k * sizeof(double) * 2);
	  tail_x_internal = (double *) blas_malloc(k * sizeof(double) * 2);
	  if (k > 0 && (head_x_internal == NULL || tail_x_internal == NULL)) {
	    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
	  }


	  FPU_FIX_START;

	  /*loop 1 */
	  xi = start_xi;
	  /* x_inti already initialized to 0 */
	  for (j = 0; j < k; j++) {

	    /* each time through loop, xi lands on next x to compute. */
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    /* preform the multiplication -
	       in this implementation we do not seperate the alpha = 1 case */
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[0] * split;
		a1 = con - x_elem[0];
		a1 = con - a1;
		a2 = x_elem[0] - a1;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		head_t1 = x_elem[0] * alpha_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[1] * split;
		a1 = con - x_elem[1];
		a1 = con - a1;
		a2 = x_elem[1] - a1;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		head_t2 = x_elem[1] * alpha_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[1] * split;
		a1 = con - x_elem[1];
		a1 = con - a1;
		a2 = x_elem[1] - a1;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		head_t1 = x_elem[1] * alpha_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[0] * split;
		a1 = con - x_elem[0];
		a1 = con - a1;
		a2 = x_elem[0] - a1;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		head_t2 = x_elem[0] * alpha_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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

	    Tij = dot_start;
	    dot_start += dot_start_inc1;

	    /*start loop buffer over in loop 1 */
	    x_inti = 0;
	    for (i = j; i > 0; i--) {
	      T_element = t_i[Tij];

	      head_temp3[0] = head_x_internal[x_inti];
	      head_temp3[1] = head_x_internal[1 + x_inti];
	      tail_temp3[0] = tail_x_internal[x_inti];
	      tail_temp3[1] = tail_x_internal[1 + x_inti];
	      {
		/* Compute complex-extra = complex-extra * real. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t, tail_t;
		head_a0 = head_temp3[0];
		tail_a0 = tail_temp3[0];
		head_a1 = head_temp3[1];
		tail_a1 = tail_temp3[1];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a0 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[0] = head_t;
		tail_temp2[0] = tail_t;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a1 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[1] = head_t;
		tail_temp2[1] = tail_t;
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
	      x_inti += inc_x_inti;
	      Tij += dot_inc;
	    }			/* for across row */


	    /* if the diagonal entry is not equal to one, then divide Xj by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element = t_i[Tij];


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
		  t1 = head_a / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_a - t12;
		  e = t11 - head_a;
		  t21 = ((-t12 - e) + (head_a - (t11 - e))) + tail_a - t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

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
		  t1 = head_a / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_a - t12;
		  e = t11 - head_a;
		  t21 = ((-t12 - e) + (head_a - (t11 - e))) + tail_a - t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_b = t1 + t2;
		  tail_b = t2 - (head_b - t1);
		}
		head_temp1[1] = head_b;
		tail_temp1[1] = tail_b;
	      }

	    }

	    /* if (diag == blas_non_unit_diag) */
	    /* place internal precision result in internal buffer */
	    head_x_internal[x_inti] = head_temp1[0];
	    tail_x_internal[x_inti] = tail_temp1[0];
	    head_x_internal[1 + x_inti] = head_temp1[1];
	    tail_x_internal[1 + x_inti] = tail_temp1[1];

	    /* place result x in same place as got x this loop */
	    x_i[xi] = head_temp1[0];
	    x_i[xi + 1] = head_temp1[1];
	    xi += incxi;
	  }			/* for j<k */
	  /*end loop 1 */


	  /* loop2 ***************************** */
	  x_inti = 0;
	  /*loop 2 continue without changing j to start */
	  for (; j < n; j++) {

	    /* each time through loop, xi lands on next x to compute. */
	    x_elem[0] = x_i[xi];
	    x_elem[1] = x_i[xi + 1];
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[0] * split;
		a1 = con - x_elem[0];
		a1 = con - a1;
		a2 = x_elem[0] - a1;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		head_t1 = x_elem[0] * alpha_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[1] * split;
		a1 = con - x_elem[1];
		a1 = con - a1;
		a2 = x_elem[1] - a1;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		head_t2 = x_elem[1] * alpha_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[1] * split;
		a1 = con - x_elem[1];
		a1 = con - a1;
		a2 = x_elem[1] - a1;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		head_t1 = x_elem[1] * alpha_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem[0] * split;
		a1 = con - x_elem[0];
		a1 = con - a1;
		a2 = x_elem[0] - a1;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		head_t2 = x_elem[0] * alpha_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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


	    Tij = dot_start;
	    dot_start += dot_start_inc2;

	    for (i = k; i > 0 && (x_inti < k_compare); i--) {
	      T_element = t_i[Tij];

	      head_temp3[0] = head_x_internal[x_inti];
	      head_temp3[1] = head_x_internal[1 + x_inti];
	      tail_temp3[0] = tail_x_internal[x_inti];
	      tail_temp3[1] = tail_x_internal[1 + x_inti];
	      {
		/* Compute complex-extra = complex-extra * real. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t, tail_t;
		head_a0 = head_temp3[0];
		tail_a0 = tail_temp3[0];
		head_a1 = head_temp3[1];
		tail_a1 = tail_temp3[1];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a0 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[0] = head_t;
		tail_temp2[0] = tail_t;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a1 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[1] = head_t;
		tail_temp2[1] = tail_t;
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
	      x_inti += inc_x_inti;
	      Tij += dot_inc;
	    }			/* for across row */
	    /*reset index to internal storage loop buffer. */
	    x_inti = 0;
	    for (; i > 0; i--) {
	      T_element = t_i[Tij];

	      head_temp3[0] = head_x_internal[x_inti];
	      head_temp3[1] = head_x_internal[1 + x_inti];
	      tail_temp3[0] = tail_x_internal[x_inti];
	      tail_temp3[1] = tail_x_internal[1 + x_inti];
	      {
		/* Compute complex-extra = complex-extra * real. */
		double head_a0, tail_a0;
		double head_a1, tail_a1;
		double head_t, tail_t;
		head_a0 = head_temp3[0];
		tail_a0 = tail_temp3[0];
		head_a1 = head_temp3[1];
		tail_a1 = tail_temp3[1];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a0 * split;
		  a11 = con - head_a0;
		  a11 = con - a11;
		  a21 = head_a0 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a0 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a0 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[0] = head_t;
		tail_temp2[0] = tail_t;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_a1 * split;
		  a11 = con - head_a1;
		  a11 = con - a11;
		  a21 = head_a1 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_a1 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_a1 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_temp2[1] = head_t;
		tail_temp2[1] = tail_t;
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
	      x_inti += inc_x_inti;
	      Tij += dot_inc;
	    }			/* for across row */


	    /* if the diagonal entry is not equal to one, then divide by 
	       the entry */
	    if (diag == blas_non_unit_diag) {
	      T_element = t_i[Tij];


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
		  t1 = head_a / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_a - t12;
		  e = t11 - head_a;
		  t21 = ((-t12 - e) + (head_a - (t11 - e))) + tail_a - t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

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
		  t1 = head_a / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_a - t12;
		  e = t11 - head_a;
		  t21 = ((-t12 - e) + (head_a - (t11 - e))) + tail_a - t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_b = t1 + t2;
		  tail_b = t2 - (head_b - t1);
		}
		head_temp1[1] = head_b;
		tail_temp1[1] = tail_b;
	      }

	    }

	    /* if (diag == blas_non_unit_diag) */
	    /* place internal precision result in internal buffer */
	    head_x_internal[x_inti] = head_temp1[0];
	    tail_x_internal[x_inti] = tail_temp1[0];
	    head_x_internal[1 + x_inti] = head_temp1[1];
	    tail_x_internal[1 + x_inti] = tail_temp1[1];
	    x_inti += inc_x_inti;
	    if (x_inti >= k_compare)
	      x_inti = 0;

	    /* place result x in same place as got x this loop */
	    x_i[xi] = head_temp1[0];
	    x_i[xi + 1] = head_temp1[1];
	    xi += incxi;
	  }			/* for j<n */

	  FPU_FIX_STOP;

	  blas_free(head_x_internal);
	  blas_free(tail_x_internal);
	}
      }
      break;
    }

  default:
    BLAS_error(routine_name, -13, prec, 0);
    break;
  }				/* end prec switch */
}				/* end BLAS_ztbsv_d_x */
