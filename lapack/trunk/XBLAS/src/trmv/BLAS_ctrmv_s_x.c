#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ctrmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const float *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec)

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
 * alpha  (input) const void*
 * 
 * T      (input) float*
 *        Triangular matrix
 *
 * ldt    (input) int 
 *        Leading dimension of T
 *
 * x      (input) const void*
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
  static const char routine_name[] = "BLAS_ctrmv_s_x";

  switch (prec) {
  case blas_prec_single:{

      int i, j;			/* used to idx matrix */
      int xj, xj0;
      int ti, tij, tij0;

      int inc_ti, inc_tij;
      int inc_x;

      const float *T_i = T;	/* internal matrix T */
      float *x_i = (float *) x;	/* internal x */
      float *alpha_i = (float *) alpha;	/* internal alpha */

      float t_elem;
      float x_elem[2];
      float prod[2];
      float sum[2];
      float tmp[2];



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



      inc_x *= 2;

      xj0 = (inc_x > 0 ? 0 : -(n - 1) * inc_x);
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	xj = xj0;
	for (j = 0; j < n; j++) {
	  x_i[xj] = 0.0;
	  x_i[xj + 1] = 0.0;
	  xj += inc_x;
	}
      } else {

	if (diag == blas_unit_diag) {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum[0] = sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < (n - 1); j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		prod[0] = x_elem[0] * t_elem;
		prod[1] = x_elem[1] * t_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    sum[0] = sum[0] + x_elem[0];
	    sum[1] = sum[1] + x_elem[1];

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj] = sum[0];
	      x_i[xj + 1] = sum[1];
	    } else {
	      {
		tmp[0] = sum[0] * alpha_i[0] - sum[1] * alpha_i[1];
		tmp[1] = sum[0] * alpha_i[1] + sum[1] * alpha_i[0];
	      }

	      x_i[xj] = tmp[0];
	      x_i[xj + 1] = tmp[1];
	    }

	    ti += inc_ti;
	  }

	} else {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum[0] = sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < n; j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		prod[0] = x_elem[0] * t_elem;
		prod[1] = x_elem[1] * t_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj - inc_x] = sum[0];
	      x_i[xj - inc_x + 1] = sum[1];
	    } else {
	      {
		tmp[0] = sum[0] * alpha_i[0] - sum[1] * alpha_i[1];
		tmp[1] = sum[0] * alpha_i[1] + sum[1] * alpha_i[0];
	      }

	      x_i[xj - inc_x] = tmp[0];
	      x_i[xj - inc_x + 1] = tmp[1];
	    }

	    ti += inc_ti;
	  }

	}

      }



      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;			/* used to idx matrix */
      int xj, xj0;
      int ti, tij, tij0;

      int inc_ti, inc_tij;
      int inc_x;

      const float *T_i = T;	/* internal matrix T */
      float *x_i = (float *) x;	/* internal x */
      float *alpha_i = (float *) alpha;	/* internal alpha */

      float t_elem;
      float x_elem[2];
      double prod[2];
      double sum[2];
      double tmp[2];



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



      inc_x *= 2;

      xj0 = (inc_x > 0 ? 0 : -(n - 1) * inc_x);
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	xj = xj0;
	for (j = 0; j < n; j++) {
	  x_i[xj] = 0.0;
	  x_i[xj + 1] = 0.0;
	  xj += inc_x;
	}
      } else {

	if (diag == blas_unit_diag) {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum[0] = sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < (n - 1); j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		prod[0] = (double) x_elem[0] * t_elem;
		prod[1] = (double) x_elem[1] * t_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    sum[0] = sum[0] + x_elem[0];
	    sum[1] = sum[1] + x_elem[1];

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj] = sum[0];
	      x_i[xj + 1] = sum[1];
	    } else {
	      {
		tmp[0] =
		  (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
		tmp[1] =
		  (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	      }
	      x_i[xj] = tmp[0];
	      x_i[xj + 1] = tmp[1];
	    }

	    ti += inc_ti;
	  }

	} else {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    sum[0] = sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < n; j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		prod[0] = (double) x_elem[0] * t_elem;
		prod[1] = (double) x_elem[1] * t_elem;
	      }
	      sum[0] = sum[0] + prod[0];
	      sum[1] = sum[1] + prod[1];

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj - inc_x] = sum[0];
	      x_i[xj - inc_x + 1] = sum[1];
	    } else {
	      {
		tmp[0] =
		  (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
		tmp[1] =
		  (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
	      }
	      x_i[xj - inc_x] = tmp[0];
	      x_i[xj - inc_x + 1] = tmp[1];
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

      const float *T_i = T;	/* internal matrix T */
      float *x_i = (float *) x;	/* internal x */
      float *alpha_i = (float *) alpha;	/* internal alpha */

      float t_elem;
      float x_elem[2];
      double head_prod[2], tail_prod[2];
      double head_sum[2], tail_sum[2];
      double head_tmp[2], tail_tmp[2];

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



      inc_x *= 2;

      xj0 = (inc_x > 0 ? 0 : -(n - 1) * inc_x);
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	xj = xj0;
	for (j = 0; j < n; j++) {
	  x_i[xj] = 0.0;
	  x_i[xj + 1] = 0.0;
	  xj += inc_x;
	}
      } else {

	if (diag == blas_unit_diag) {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < (n - 1); j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		head_prod[0] = (double) x_elem[0] * t_elem;
		tail_prod[0] = 0.0;
		head_prod[1] = (double) x_elem[1] * t_elem;
		tail_prod[1] = 0.0;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    {
	      double cd[2];
	      cd[0] = (double) x_elem[0];
	      cd[1] = (double) x_elem[1];
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		{
		  /* Compute double-double = double-double + double. */
		  double e, t1, t2;

		  /* Knuth trick. */
		  t1 = head_a + cd[0];
		  e = t1 - head_a;
		  t2 = ((cd[0] - e) + (head_a - (t1 - e))) + tail_a;

		  /* The result is t1 + t2, after normalization. */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		{
		  /* Compute double-double = double-double + double. */
		  double e, t1, t2;

		  /* Knuth trick. */
		  t1 = head_a + cd[1];
		  e = t1 - head_a;
		  t2 = ((cd[1] - e) + (head_a - (t1 - e))) + tail_a;

		  /* The result is t1 + t2, after normalization. */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }
	    }

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj] = head_sum[0];
	      x_i[xj + 1] = head_sum[1];
	    } else {
	      {
		double cd[2];
		cd[0] = (double) alpha_i[0];
		cd[1] = (double) alpha_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
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
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
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
		  head_tmp[0] = head_t1;
		  tail_tmp[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
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
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
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
		  head_tmp[1] = head_t1;
		  tail_tmp[1] = tail_t1;
		}

	      }
	      x_i[xj] = head_tmp[0];
	      x_i[xj + 1] = head_tmp[1];
	    }

	    ti += inc_ti;
	  }

	} else {


	  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	  tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	  for (i = 0; i < n; i++) {

	    head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;

	    xj = xj0;
	    tij = ti + tij0;
	    for (j = i; j < n; j++) {

	      t_elem = T_i[tij];

	      x_elem[0] = x_i[xj];
	      x_elem[1] = x_i[xj + 1];
	      {
		head_prod[0] = (double) x_elem[0] * t_elem;
		tail_prod[0] = 0.0;
		head_prod[1] = (double) x_elem[1] * t_elem;
		tail_prod[1] = 0.0;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sum[0];
		tail_a = tail_sum[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[0] = head_t;
		tail_sum[0] = tail_t;
		/* Imaginary part */
		head_a = head_sum[1];
		tail_a = tail_sum[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sum[1] = head_t;
		tail_sum[1] = tail_t;
	      }

	      xj += inc_x;
	      tij += inc_tij;
	    }

	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      x_i[xj - inc_x] = head_sum[0];
	      x_i[xj - inc_x + 1] = head_sum[1];
	    } else {
	      {
		double cd[2];
		cd[0] = (double) alpha_i[0];
		cd[1] = (double) alpha_i[1];
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
		  /* real part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a0 * split;
		    a11 = con - head_a0;
		    a11 = con - a11;
		    a21 = head_a0 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a0 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[0];
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
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a1 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[1];
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
		  head_tmp[0] = head_t1;
		  tail_tmp[0] = tail_t1;
		  /* imaginary part */
		  {
		    /* Compute double-double = double-double * double. */
		    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		    con = head_a1 * split;
		    a11 = con - head_a1;
		    a11 = con - a11;
		    a21 = head_a1 - a11;
		    con = cd[0] * split;
		    b1 = con - cd[0];
		    b1 = con - b1;
		    b2 = cd[0] - b1;

		    c11 = head_a1 * cd[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a1 * cd[0];
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
		    con = cd[1] * split;
		    b1 = con - cd[1];
		    b1 = con - b1;
		    b2 = cd[1] - b1;

		    c11 = head_a0 * cd[1];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		    c2 = tail_a0 * cd[1];
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
		  head_tmp[1] = head_t1;
		  tail_tmp[1] = tail_t1;
		}

	      }
	      x_i[xj - inc_x] = head_tmp[0];
	      x_i[xj - inc_x + 1] = head_tmp[1];
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
