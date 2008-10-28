#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_dtrmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  double alpha, const double *T, int ldt,
		  double *x, int incx, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes x <-- alpha * T * x, where T is a triangular matrix.
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
 * alpha  (input) double
 * 
 * T      (input) double*
 *        Triangular matrix
 *
 * ldt    (input) int 
 *        Leading dimension of T
 *
 * x      (input) const double*
 *    Array of length n.
 * 
 * incx   (input) int
 *     The stride used to access components x[i].
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
  static const char routine_name[] = "BLAS_dtrmv_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;			/* used to idx matrix */
      int xj, xj0;
      int ti, tij, tij0;

      int inc_ti, inc_tij;
      int inc_x;

      const double *T_i = T;	/* internal matrix T */
      double *x_i = x;		/* internal x */
      double alpha_i = alpha;	/* internal alpha */

      double t_elem;
      double x_elem;
      double prod;
      double sum;
      double tmp;



      /* all error calls */
      if ((order != blas_rowmajor && order != blas_colmajor) ||
	  (uplo != blas_upper && uplo != blas_lower) ||
	  (trans != blas_trans &&
	   trans != blas_no_trans &&
	   trans != blas_conj_trans) ||
	  (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
	  (ldt < n) || (incx == 0)) {
	BLAS_error(routine_name, 0, 0, NULL);
      } else if (n <= 0) {
	BLAS_error(routine_name, -4, n, NULL);
      } else if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }

      if (trans == blas_no_trans) {
	if (uplo == blas_upper) {
	  inc_x = -incx;
	  if (order == blas_rowmajor) {
	    inc_ti = ldt;
	    inc_tij = -1;
	  } else {
	    inc_ti = 1;
	    inc_tij = -ldt;
	  }
	} else {
	  inc_x = incx;
	  if (order == blas_rowmajor) {
	    inc_ti = -ldt;
	    inc_tij = 1;
	  } else {
	    inc_ti = -1;
	    inc_tij = ldt;
	  }
	}
      } else {
	if (uplo == blas_upper) {
	  inc_x = incx;
	  if (order == blas_rowmajor) {
	    inc_ti = -1;
	    inc_tij = ldt;
	  } else {
	    inc_ti = -ldt;
	    inc_tij = 1;
	  }
	} else {
	  inc_x = -incx;
	  if (order == blas_rowmajor) {
	    inc_ti = 1;
	    inc_tij = -ldt;
	  } else {
	    inc_ti = ldt;
	    inc_tij = -1;
	  }
	}
      }





      xj0 = (inc_x > 0 ? 0 : -(n - 1) * inc_x);
      if (alpha_i == 0.0) {
	xj = xj0;
	for (j = 0; j < n; j++) {
	  x_i[xj] = 0.0;
	  xj += inc_x;
	}
      } else {

	if (diag == blas_unit_diag) {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < (n - 1); j++) {

	      t_elem = T_i[tij];

	      x_elem = x_i[xj];
	      prod = x_elem * t_elem;
	      sum = sum + prod;

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    x_elem = x_i[xj];
	    sum = sum + x_elem;

	    if (alpha_i == 1.0) {
	      x_i[xj] = sum;
	    } else {
	      tmp = sum * alpha_i;
	      x_i[xj] = tmp;
	    }

	    ti += inc_ti;
	  }

	} else {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < n; j++) {

	      t_elem = T_i[tij];

	      x_elem = x_i[xj];
	      prod = x_elem * t_elem;
	      sum = sum + prod;

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    if (alpha_i == 1.0) {
	      x_i[xj - inc_x] = sum;
	    } else {
	      tmp = sum * alpha_i;
	      x_i[xj - inc_x] = tmp;
	    }

	    ti += inc_ti;
	  }

	}

      }



      break;
    }

  case blas_prec_extra:{

      int i, j;			/* used to idx matrix */
      int xj, xj0;
      int ti, tij, tij0;

      int inc_ti, inc_tij;
      int inc_x;

      const double *T_i = T;	/* internal matrix T */
      double *x_i = x;		/* internal x */
      double alpha_i = alpha;	/* internal alpha */

      double t_elem;
      double x_elem;
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp, tail_tmp;

      FPU_FIX_DECL;
      FPU_FIX_START;

      /* all error calls */
      if ((order != blas_rowmajor && order != blas_colmajor) ||
	  (uplo != blas_upper && uplo != blas_lower) ||
	  (trans != blas_trans &&
	   trans != blas_no_trans &&
	   trans != blas_conj_trans) ||
	  (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
	  (ldt < n) || (incx == 0)) {
	BLAS_error(routine_name, 0, 0, NULL);
      } else if (n <= 0) {
	BLAS_error(routine_name, -4, n, NULL);
      } else if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }

      if (trans == blas_no_trans) {
	if (uplo == blas_upper) {
	  inc_x = -incx;
	  if (order == blas_rowmajor) {
	    inc_ti = ldt;
	    inc_tij = -1;
	  } else {
	    inc_ti = 1;
	    inc_tij = -ldt;
	  }
	} else {
	  inc_x = incx;
	  if (order == blas_rowmajor) {
	    inc_ti = -ldt;
	    inc_tij = 1;
	  } else {
	    inc_ti = -1;
	    inc_tij = ldt;
	  }
	}
      } else {
	if (uplo == blas_upper) {
	  inc_x = incx;
	  if (order == blas_rowmajor) {
	    inc_ti = -1;
	    inc_tij = ldt;
	  } else {
	    inc_ti = -ldt;
	    inc_tij = 1;
	  }
	} else {
	  inc_x = -incx;
	  if (order == blas_rowmajor) {
	    inc_ti = 1;
	    inc_tij = -ldt;
	  } else {
	    inc_ti = ldt;
	    inc_tij = -1;
	  }
	}
      }





      xj0 = (inc_x > 0 ? 0 : -(n - 1) * inc_x);
      if (alpha_i == 0.0) {
	xj = xj0;
	for (j = 0; j < n; j++) {
	  x_i[xj] = 0.0;
	  xj += inc_x;
	}
      } else {

	if (diag == blas_unit_diag) {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    head_sum = tail_sum = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < (n - 1); j++) {

	      t_elem = T_i[tij];

	      x_elem = x_i[xj];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = t_elem * split;
		b1 = con - t_elem;
		b1 = con - b1;
		b2 = t_elem - b1;

		head_prod = x_elem * t_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    x_elem = x_i[xj];
	    {
	      /* Compute double-double = double-double + double. */
	      double e, t1, t2;

	      /* Knuth trick. */
	      t1 = head_sum + x_elem;
	      e = t1 - head_sum;
	      t2 = ((x_elem - e) + (head_sum - (t1 - e))) + tail_sum;

	      /* The result is t1 + t2, after normalization. */
	      head_sum = t1 + t2;
	      tail_sum = t2 - (head_sum - t1);
	    }

	    if (alpha_i == 1.0) {
	      x_i[xj] = head_sum;
	    } else {
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_sum * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_tmp = t1 + t2;
		tail_tmp = t2 - (head_tmp - t1);
	      }
	      x_i[xj] = head_tmp;
	    }

	    ti += inc_ti;
	  }

	} else {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    head_sum = tail_sum = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < n; j++) {

	      t_elem = T_i[tij];

	      x_elem = x_i[xj];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = x_elem * split;
		a1 = con - x_elem;
		a1 = con - a1;
		a2 = x_elem - a1;
		con = t_elem * split;
		b1 = con - t_elem;
		b1 = con - b1;
		b2 = t_elem - b1;

		head_prod = x_elem * t_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    if (alpha_i == 1.0) {
	      x_i[xj - inc_x] = head_sum;
	    } else {
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_sum * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_tmp = t1 + t2;
		tail_tmp = t2 - (head_tmp - t1);
	      }
	      x_i[xj - inc_x] = head_tmp;
	    }

	    ti += inc_ti;
	  }

	}

      }

      FPU_FIX_STOP;

      break;
    }
  }
}
