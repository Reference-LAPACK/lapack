#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ztrmv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const void *T, int ldt,
		  void *x, int incx)

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
 * T      (input) void*
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
 */
{
  static const char routine_name[] = "BLAS_ztrmv_c";

  int i, j;			/* used to idx matrix */
  int xj, xj0;
  int ti, tij, tij0;

  int inc_ti, inc_tij;
  int inc_x;

  const float *T_i = (float *) T;	/* internal matrix T */
  double *x_i = (double *) x;	/* internal x */
  double *alpha_i = (double *) alpha;	/* internal alpha */

  float t_elem[2];
  double x_elem[2];
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

  inc_ti *= 2;
  inc_tij *= 2;
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
      if (trans == blas_conj_trans) {


	ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	for (i = 0; i < n; i++) {

	  sum[0] = sum[1] = 0.0;

	  xj = xj0;
	  tij = ti + tij0;
	  for (j = i; j < (n - 1); j++) {

	    t_elem[0] = T_i[tij];
	    t_elem[1] = T_i[tij + 1];
	    t_elem[1] = -t_elem[1];
	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    {
	      prod[0] =
		(double) x_elem[0] * t_elem[0] -
		(double) x_elem[1] * t_elem[1];
	      prod[1] =
		(double) x_elem[0] * t_elem[1] +
		(double) x_elem[1] * t_elem[0];
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
	  for (j = i; j < (n - 1); j++) {

	    t_elem[0] = T_i[tij];
	    t_elem[1] = T_i[tij + 1];

	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    {
	      prod[0] =
		(double) x_elem[0] * t_elem[0] -
		(double) x_elem[1] * t_elem[1];
	      prod[1] =
		(double) x_elem[0] * t_elem[1] +
		(double) x_elem[1] * t_elem[0];
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

      }
    } else {
      if (trans == blas_conj_trans) {


	ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	for (i = 0; i < n; i++) {

	  sum[0] = sum[1] = 0.0;

	  xj = xj0;
	  tij = ti + tij0;
	  for (j = i; j < n; j++) {

	    t_elem[0] = T_i[tij];
	    t_elem[1] = T_i[tij + 1];
	    t_elem[1] = -t_elem[1];
	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    {
	      prod[0] =
		(double) x_elem[0] * t_elem[0] -
		(double) x_elem[1] * t_elem[1];
	      prod[1] =
		(double) x_elem[0] * t_elem[1] +
		(double) x_elem[1] * t_elem[0];
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

      } else {


	ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
	tij0 = (inc_tij > 0 ? 0 : -(n - 1) * inc_tij);
	for (i = 0; i < n; i++) {

	  sum[0] = sum[1] = 0.0;

	  xj = xj0;
	  tij = ti + tij0;
	  for (j = i; j < n; j++) {

	    t_elem[0] = T_i[tij];
	    t_elem[1] = T_i[tij + 1];

	    x_elem[0] = x_i[xj];
	    x_elem[1] = x_i[xj + 1];
	    {
	      prod[0] =
		(double) x_elem[0] * t_elem[0] -
		(double) x_elem[1] * t_elem[1];
	      prod[1] =
		(double) x_elem[0] * t_elem[1] +
		(double) x_elem[1] * t_elem[0];
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

  }



}
