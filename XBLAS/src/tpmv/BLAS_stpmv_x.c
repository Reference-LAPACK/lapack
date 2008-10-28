#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_stpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, float alpha, const float *tp,
		  float *x, int incx, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes x = alpha * tp * x, x = alpha * tp_transpose * x,
 * or x = alpha * tp_conjugate_transpose where tp is a triangular
 * packed matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *
 * diag         (input) blas_diag_type
 *              Whether the diagonal entries of tp are 1
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input) float
 *
 * tp           (input) float*
 *
 * x            (input) float*
 *
 * incx         (input) int
 *              The stride for vector x.
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
  static const char routine_name[] = "BLAS_stpmv_x";


  switch (prec) {
  case blas_prec_single:{

      {
	int matrix_row, step, tp_index, tp_start, x_index, x_start;
	int inctp, x_index2, stride, col_index, inctp2;

	float alpha_i = alpha;

	const float *tp_i = tp;
	float *x_i = x;
	float rowsum;
	float rowtmp;
	float result;
	float matval;
	float vecval;
	float one;


	one = 1.0;

	inctp = 1;



	if (incx < 0)
	  x_start = (-n + 1) * incx;
	else
	  x_start = 0;

	if (n < 1) {
	  return;
	}

	/* Check for error conditions. */
	if (order != blas_colmajor && order != blas_rowmajor) {
	  BLAS_error(routine_name, -1, order, NULL);
	}
	if (uplo != blas_upper && uplo != blas_lower) {
	  BLAS_error(routine_name, -2, uplo, NULL);
	}
	if (incx == 0) {
	  BLAS_error(routine_name, -9, incx, NULL);
	}



	{
	  if ((uplo == blas_upper &&
	       trans == blas_no_trans && order == blas_rowmajor) ||
	      (uplo == blas_lower &&
	       trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = 0;
	    tp_index = tp_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      x_index = x_start + incx * matrix_row;
	      x_index2 = x_index;
	      col_index = matrix_row;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      while (col_index < n) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  rowtmp = vecval * one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = matval * vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index2] = result;
	    }
	  } else if ((uplo == blas_upper &&
		      trans == blas_no_trans && order == blas_colmajor) ||
		     (uplo == blas_lower &&
		      trans != blas_no_trans && order == blas_rowmajor)) {
	    tp_start = ((n - 1) * n) / 2;
	    inctp2 = n - 1;
	    x_index2 = x_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++, inctp2 = n - 1) {
	      x_index = x_start + incx * (n - 1);
	      tp_index = (tp_start + matrix_row) * inctp;
	      col_index = (n - 1) - matrix_row;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      while (col_index >= 0) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  rowtmp = vecval * one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = matval * vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index2] = result;
	      x_index2 += incx;
	    }
	  } else if ((uplo == blas_lower &&
		      trans == blas_no_trans && order == blas_rowmajor) ||
		     (uplo == blas_upper &&
		      trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = (n - 1) + ((n - 1) * n) / 2;
	    tp_index = tp_start * inctp;
	    x_index = x_start + (n - 1) * incx;

	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      x_index2 = x_index;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  rowtmp = vecval * one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = matval * vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index2 -= incx;
		tp_index -= inctp;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index] = result;
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  rowtmp = vecval * one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = matval * vecval;
		}
		rowsum = rowsum + rowtmp;
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index] = result;
	      x_index -= incx;
	    }
	  }
	}


      }
      break;
    }
  case blas_prec_double:
  case blas_prec_indigenous:{

      {
	int matrix_row, step, tp_index, tp_start, x_index, x_start;
	int inctp, x_index2, stride, col_index, inctp2;

	float alpha_i = alpha;

	const float *tp_i = tp;
	float *x_i = x;
	double rowsum;
	double rowtmp;
	double result;
	float matval;
	float vecval;
	float one;


	one = 1.0;

	inctp = 1;



	if (incx < 0)
	  x_start = (-n + 1) * incx;
	else
	  x_start = 0;

	if (n < 1) {
	  return;
	}

	/* Check for error conditions. */
	if (order != blas_colmajor && order != blas_rowmajor) {
	  BLAS_error(routine_name, -1, order, NULL);
	}
	if (uplo != blas_upper && uplo != blas_lower) {
	  BLAS_error(routine_name, -2, uplo, NULL);
	}
	if (incx == 0) {
	  BLAS_error(routine_name, -9, incx, NULL);
	}



	{
	  if ((uplo == blas_upper &&
	       trans == blas_no_trans && order == blas_rowmajor) ||
	      (uplo == blas_lower &&
	       trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = 0;
	    tp_index = tp_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      x_index = x_start + incx * matrix_row;
	      x_index2 = x_index;
	      col_index = matrix_row;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      while (col_index < n) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  rowtmp = (double) vecval *one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = (double) matval *vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index2] = result;
	    }
	  } else if ((uplo == blas_upper &&
		      trans == blas_no_trans && order == blas_colmajor) ||
		     (uplo == blas_lower &&
		      trans != blas_no_trans && order == blas_rowmajor)) {
	    tp_start = ((n - 1) * n) / 2;
	    inctp2 = n - 1;
	    x_index2 = x_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++, inctp2 = n - 1) {
	      x_index = x_start + incx * (n - 1);
	      tp_index = (tp_start + matrix_row) * inctp;
	      col_index = (n - 1) - matrix_row;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      while (col_index >= 0) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  rowtmp = (double) vecval *one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = (double) matval *vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index2] = result;
	      x_index2 += incx;
	    }
	  } else if ((uplo == blas_lower &&
		      trans == blas_no_trans && order == blas_rowmajor) ||
		     (uplo == blas_upper &&
		      trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = (n - 1) + ((n - 1) * n) / 2;
	    tp_index = tp_start * inctp;
	    x_index = x_start + (n - 1) * incx;

	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      x_index2 = x_index;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  rowtmp = (double) vecval *one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = (double) matval *vecval;
		}
		rowsum = rowsum + rowtmp;
		x_index2 -= incx;
		tp_index -= inctp;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index] = result;
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      rowsum = 0.0;
	      rowtmp = 0.0;
	      result = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  rowtmp = (double) vecval *one;
		} else {
		  matval = tp_i[tp_index];
		  rowtmp = (double) matval *vecval;
		}
		rowsum = rowsum + rowtmp;
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
	      result = rowsum * alpha_i;
	      x_i[x_index] = result;
	      x_index -= incx;
	    }
	  }
	}


      }
      break;
    }

  case blas_prec_extra:{

      {
	int matrix_row, step, tp_index, tp_start, x_index, x_start;
	int inctp, x_index2, stride, col_index, inctp2;

	float alpha_i = alpha;

	const float *tp_i = tp;
	float *x_i = x;
	double head_rowsum, tail_rowsum;
	double head_rowtmp, tail_rowtmp;
	double head_result, tail_result;
	float matval;
	float vecval;
	float one;

	FPU_FIX_DECL;
	one = 1.0;

	inctp = 1;



	if (incx < 0)
	  x_start = (-n + 1) * incx;
	else
	  x_start = 0;

	if (n < 1) {
	  return;
	}

	/* Check for error conditions. */
	if (order != blas_colmajor && order != blas_rowmajor) {
	  BLAS_error(routine_name, -1, order, NULL);
	}
	if (uplo != blas_upper && uplo != blas_lower) {
	  BLAS_error(routine_name, -2, uplo, NULL);
	}
	if (incx == 0) {
	  BLAS_error(routine_name, -9, incx, NULL);
	}
	FPU_FIX_START;


	{
	  if ((uplo == blas_upper &&
	       trans == blas_no_trans && order == blas_rowmajor) ||
	      (uplo == blas_lower &&
	       trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = 0;
	    tp_index = tp_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      x_index = x_start + incx * matrix_row;
	      x_index2 = x_index;
	      col_index = matrix_row;
	      head_rowsum = tail_rowsum = 0.0;
	      head_rowtmp = tail_rowtmp = 0.0;
	      head_result = tail_result = 0.0;
	      while (col_index < n) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  head_rowtmp = (double) vecval *one;
		  tail_rowtmp = 0.0;
		} else {
		  matval = tp_i[tp_index];
		  head_rowtmp = (double) matval *vecval;
		  tail_rowtmp = 0.0;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_rowsum + head_rowtmp;
		  bv = s1 - head_rowsum;
		  s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_rowsum + tail_rowtmp;
		  bv = t1 - tail_rowsum;
		  t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_rowsum = t1 + t2;
		  tail_rowsum = t2 - (head_rowsum - t1);
		}
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
	      {
		double dt = (double) alpha_i;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_rowsum * split;
		  a11 = con - head_rowsum;
		  a11 = con - a11;
		  a21 = head_rowsum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_rowsum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_rowsum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_result = t1 + t2;
		  tail_result = t2 - (head_result - t1);
		}
	      }
	      x_i[x_index2] = head_result;
	    }
	  } else if ((uplo == blas_upper &&
		      trans == blas_no_trans && order == blas_colmajor) ||
		     (uplo == blas_lower &&
		      trans != blas_no_trans && order == blas_rowmajor)) {
	    tp_start = ((n - 1) * n) / 2;
	    inctp2 = n - 1;
	    x_index2 = x_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++, inctp2 = n - 1) {
	      x_index = x_start + incx * (n - 1);
	      tp_index = (tp_start + matrix_row) * inctp;
	      col_index = (n - 1) - matrix_row;
	      head_rowsum = tail_rowsum = 0.0;
	      head_rowtmp = tail_rowtmp = 0.0;
	      head_result = tail_result = 0.0;
	      while (col_index >= 0) {
		vecval = x_i[x_index];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  head_rowtmp = (double) vecval *one;
		  tail_rowtmp = 0.0;
		} else {
		  matval = tp_i[tp_index];
		  head_rowtmp = (double) matval *vecval;
		  tail_rowtmp = 0.0;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_rowsum + head_rowtmp;
		  bv = s1 - head_rowsum;
		  s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_rowsum + tail_rowtmp;
		  bv = t1 - tail_rowsum;
		  t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_rowsum = t1 + t2;
		  tail_rowsum = t2 - (head_rowsum - t1);
		}
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
	      {
		double dt = (double) alpha_i;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_rowsum * split;
		  a11 = con - head_rowsum;
		  a11 = con - a11;
		  a21 = head_rowsum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_rowsum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_rowsum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_result = t1 + t2;
		  tail_result = t2 - (head_result - t1);
		}
	      }
	      x_i[x_index2] = head_result;
	      x_index2 += incx;
	    }
	  } else if ((uplo == blas_lower &&
		      trans == blas_no_trans && order == blas_rowmajor) ||
		     (uplo == blas_upper &&
		      trans != blas_no_trans && order == blas_colmajor)) {
	    tp_start = (n - 1) + ((n - 1) * n) / 2;
	    tp_index = tp_start * inctp;
	    x_index = x_start + (n - 1) * incx;

	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      x_index2 = x_index;
	      head_rowsum = tail_rowsum = 0.0;
	      head_rowtmp = tail_rowtmp = 0.0;
	      head_result = tail_result = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  head_rowtmp = (double) vecval *one;
		  tail_rowtmp = 0.0;
		} else {
		  matval = tp_i[tp_index];
		  head_rowtmp = (double) matval *vecval;
		  tail_rowtmp = 0.0;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_rowsum + head_rowtmp;
		  bv = s1 - head_rowsum;
		  s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_rowsum + tail_rowtmp;
		  bv = t1 - tail_rowsum;
		  t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_rowsum = t1 + t2;
		  tail_rowsum = t2 - (head_rowsum - t1);
		}
		x_index2 -= incx;
		tp_index -= inctp;
	      }
	      {
		double dt = (double) alpha_i;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_rowsum * split;
		  a11 = con - head_rowsum;
		  a11 = con - a11;
		  a21 = head_rowsum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_rowsum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_rowsum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_result = t1 + t2;
		  tail_result = t2 - (head_result - t1);
		}
	      }
	      x_i[x_index] = head_result;
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      head_rowsum = tail_rowsum = 0.0;
	      head_rowtmp = tail_rowtmp = 0.0;
	      head_result = tail_result = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval = x_i[x_index2];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  head_rowtmp = (double) vecval *one;
		  tail_rowtmp = 0.0;
		} else {
		  matval = tp_i[tp_index];
		  head_rowtmp = (double) matval *vecval;
		  tail_rowtmp = 0.0;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_rowsum + head_rowtmp;
		  bv = s1 - head_rowsum;
		  s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_rowsum + tail_rowtmp;
		  bv = t1 - tail_rowsum;
		  t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_rowsum = t1 + t2;
		  tail_rowsum = t2 - (head_rowsum - t1);
		}
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
	      {
		double dt = (double) alpha_i;
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_rowsum * split;
		  a11 = con - head_rowsum;
		  a11 = con - a11;
		  a21 = head_rowsum - a11;
		  con = dt * split;
		  b1 = con - dt;
		  b1 = con - b1;
		  b2 = dt - b1;

		  c11 = head_rowsum * dt;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_rowsum * dt;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_result = t1 + t2;
		  tail_result = t2 - (head_result - t1);
		}
	      }
	      x_i[x_index] = head_result;
	      x_index -= incx;
	    }
	  }
	}


      }
      break;
    }
  }

}
