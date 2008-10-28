#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ctpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp,
		  void *x, int incx, enum blas_prec_type prec)

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
 * alpha        (input) const void*
 *
 * tp           (input) void*
 *
 * x            (input) void*
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
  static const char routine_name[] = "BLAS_ctpmv_x";


  switch (prec) {
  case blas_prec_single:{

      {
	int matrix_row, step, tp_index, tp_start, x_index, x_start;
	int inctp, x_index2, stride, col_index, inctp2;

	float *alpha_i = (float *) alpha;

	const float *tp_i = (float *) tp;
	float *x_i = (float *) x;
	float rowsum[2];
	float rowtmp[2];
	float result[2];
	float matval[2];
	float vecval[2];
	float one[2];


	one[0] = 1.0;
	one[1] = 0.0;

	inctp = 1;
	inctp *= 2;
	incx *= 2;

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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      while (col_index < n) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  {
		    rowtmp[0] = vecval[0] * one[0] - vecval[1] * one[1];
		    rowtmp[1] = vecval[0] * one[1] + vecval[1] * one[0];
		  }

		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] = matval[0] * vecval[0] - matval[1] * vecval[1];
		    rowtmp[1] = matval[0] * vecval[1] + matval[1] * vecval[0];
		  }

		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
	      {
		result[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		result[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
	      }

	      x_i[x_index2] = result[0];
	      x_i[x_index2 + 1] = result[1];
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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      while (col_index >= 0) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  {
		    rowtmp[0] = vecval[0] * one[0] - vecval[1] * one[1];
		    rowtmp[1] = vecval[0] * one[1] + vecval[1] * one[0];
		  }

		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] = matval[0] * vecval[0] - matval[1] * vecval[1];
		    rowtmp[1] = matval[0] * vecval[1] + matval[1] * vecval[0];
		  }

		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
	      {
		result[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		result[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
	      }

	      x_i[x_index2] = result[0];
	      x_i[x_index2 + 1] = result[1];
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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  {
		    rowtmp[0] = vecval[0] * one[0] - vecval[1] * one[1];
		    rowtmp[1] = vecval[0] * one[1] + vecval[1] * one[0];
		  }

		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] = matval[0] * vecval[0] - matval[1] * vecval[1];
		    rowtmp[1] = matval[0] * vecval[1] + matval[1] * vecval[0];
		  }

		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index2 -= incx;
		tp_index -= inctp;
	      }
	      {
		result[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		result[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
	      }

	      x_i[x_index] = result[0];
	      x_i[x_index + 1] = result[1];
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  {
		    rowtmp[0] = vecval[0] * one[0] - vecval[1] * one[1];
		    rowtmp[1] = vecval[0] * one[1] + vecval[1] * one[0];
		  }

		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] = matval[0] * vecval[0] - matval[1] * vecval[1];
		    rowtmp[1] = matval[0] * vecval[1] + matval[1] * vecval[0];
		  }

		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
	      {
		result[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		result[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
	      }

	      x_i[x_index] = result[0];
	      x_i[x_index + 1] = result[1];
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

	float *alpha_i = (float *) alpha;

	const float *tp_i = (float *) tp;
	float *x_i = (float *) x;
	double rowsum[2];
	double rowtmp[2];
	double result[2];
	float matval[2];
	float vecval[2];
	float one[2];


	one[0] = 1.0;
	one[1] = 0.0;

	inctp = 1;
	inctp *= 2;
	incx *= 2;

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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      while (col_index < n) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  {
		    rowtmp[0] =
		      (double) vecval[0] * one[0] -
		      (double) vecval[1] * one[1];
		    rowtmp[1] =
		      (double) vecval[0] * one[1] +
		      (double) vecval[1] * one[0];
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
	      {
		result[0] =
		  (double) rowsum[0] * alpha_i[0] -
		  (double) rowsum[1] * alpha_i[1];
		result[1] =
		  (double) rowsum[0] * alpha_i[1] +
		  (double) rowsum[1] * alpha_i[0];
	      }
	      x_i[x_index2] = result[0];
	      x_i[x_index2 + 1] = result[1];
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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      while (col_index >= 0) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  {
		    rowtmp[0] =
		      (double) vecval[0] * one[0] -
		      (double) vecval[1] * one[1];
		    rowtmp[1] =
		      (double) vecval[0] * one[1] +
		      (double) vecval[1] * one[0];
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
	      {
		result[0] =
		  (double) rowsum[0] * alpha_i[0] -
		  (double) rowsum[1] * alpha_i[1];
		result[1] =
		  (double) rowsum[0] * alpha_i[1] +
		  (double) rowsum[1] * alpha_i[0];
	      }
	      x_i[x_index2] = result[0];
	      x_i[x_index2 + 1] = result[1];
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
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  {
		    rowtmp[0] =
		      (double) vecval[0] * one[0] -
		      (double) vecval[1] * one[1];
		    rowtmp[1] =
		      (double) vecval[0] * one[1] +
		      (double) vecval[1] * one[0];
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		x_index2 -= incx;
		tp_index -= inctp;
	      }
	      {
		result[0] =
		  (double) rowsum[0] * alpha_i[0] -
		  (double) rowsum[1] * alpha_i[1];
		result[1] =
		  (double) rowsum[0] * alpha_i[1] +
		  (double) rowsum[1] * alpha_i[0];
	      }
	      x_i[x_index] = result[0];
	      x_i[x_index + 1] = result[1];
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      rowsum[0] = rowsum[1] = 0.0;
	      rowtmp[0] = rowtmp[1] = 0.0;
	      result[0] = result[1] = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  {
		    rowtmp[0] =
		      (double) vecval[0] * one[0] -
		      (double) vecval[1] * one[1];
		    rowtmp[1] =
		      (double) vecval[0] * one[1] +
		      (double) vecval[1] * one[0];
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		}
		rowsum[0] = rowsum[0] + rowtmp[0];
		rowsum[1] = rowsum[1] + rowtmp[1];
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
	      {
		result[0] =
		  (double) rowsum[0] * alpha_i[0] -
		  (double) rowsum[1] * alpha_i[1];
		result[1] =
		  (double) rowsum[0] * alpha_i[1] +
		  (double) rowsum[1] * alpha_i[0];
	      }
	      x_i[x_index] = result[0];
	      x_i[x_index + 1] = result[1];
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

	float *alpha_i = (float *) alpha;

	const float *tp_i = (float *) tp;
	float *x_i = (float *) x;
	double head_rowsum[2], tail_rowsum[2];
	double head_rowtmp[2], tail_rowtmp[2];
	double head_result[2], tail_result[2];
	float matval[2];
	float vecval[2];
	float one[2];

	FPU_FIX_DECL;
	one[0] = 1.0;
	one[1] = 0.0;

	inctp = 1;
	inctp *= 2;
	incx *= 2;

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
	      head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		tail_rowsum[1] = 0.0;
	      head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		tail_rowtmp[1] = 0.0;
	      head_result[0] = head_result[1] = tail_result[0] =
		tail_result[1] = 0.0;
	      while (col_index < n) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) vecval[0] * one[0];
		    d2 = (double) -vecval[1] * one[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) vecval[0] * one[1];
		    d2 = (double) vecval[1] * one[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) matval[0] * vecval[0];
		    d2 = (double) -matval[1] * vecval[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) matval[0] * vecval[1];
		    d2 = (double) matval[1] * vecval[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_rowsum[0];
		  tail_a = tail_rowsum[0];
		  head_b = head_rowtmp[0];
		  tail_b = tail_rowtmp[0];
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
		  head_rowsum[0] = head_t;
		  tail_rowsum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_rowsum[1];
		  tail_a = tail_rowsum[1];
		  head_b = head_rowtmp[1];
		  tail_b = tail_rowtmp[1];
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
		  head_rowsum[1] = head_t;
		  tail_rowsum[1] = tail_t;
		}
		x_index += incx;
		tp_index += inctp;
		col_index++;
	      }
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
		  head_a0 = head_rowsum[0];
		  tail_a0 = tail_rowsum[0];
		  head_a1 = head_rowsum[1];
		  tail_a1 = tail_rowsum[1];
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
		  head_result[0] = head_t1;
		  tail_result[0] = tail_t1;
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
		  head_result[1] = head_t1;
		  tail_result[1] = tail_t1;
		}

	      }
	      x_i[x_index2] = head_result[0];
	      x_i[x_index2 + 1] = head_result[1];
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
	      head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		tail_rowsum[1] = 0.0;
	      head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		tail_rowtmp[1] = 0.0;
	      head_result[0] = head_result[1] = tail_result[0] =
		tail_result[1] = 0.0;
	      while (col_index >= 0) {
		vecval[0] = x_i[x_index];
		vecval[1] = x_i[x_index + 1];
		if ((diag == blas_unit_diag) && (col_index == 0)) {
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) vecval[0] * one[0];
		    d2 = (double) -vecval[1] * one[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) vecval[0] * one[1];
		    d2 = (double) vecval[1] * one[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) matval[0] * vecval[0];
		    d2 = (double) -matval[1] * vecval[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) matval[0] * vecval[1];
		    d2 = (double) matval[1] * vecval[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_rowsum[0];
		  tail_a = tail_rowsum[0];
		  head_b = head_rowtmp[0];
		  tail_b = tail_rowtmp[0];
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
		  head_rowsum[0] = head_t;
		  tail_rowsum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_rowsum[1];
		  tail_a = tail_rowsum[1];
		  head_b = head_rowtmp[1];
		  tail_b = tail_rowtmp[1];
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
		  head_rowsum[1] = head_t;
		  tail_rowsum[1] = tail_t;
		}
		x_index -= incx;
		tp_index -= inctp2 * inctp;
		inctp2--;
		col_index--;
	      }
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
		  head_a0 = head_rowsum[0];
		  tail_a0 = tail_rowsum[0];
		  head_a1 = head_rowsum[1];
		  tail_a1 = tail_rowsum[1];
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
		  head_result[0] = head_t1;
		  tail_result[0] = tail_t1;
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
		  head_result[1] = head_t1;
		  tail_result[1] = tail_t1;
		}

	      }
	      x_i[x_index2] = head_result[0];
	      x_i[x_index2 + 1] = head_result[1];
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
	      head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		tail_rowsum[1] = 0.0;
	      head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		tail_rowtmp[1] = 0.0;
	      head_result[0] = head_result[1] = tail_result[0] =
		tail_result[1] = 0.0;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == 0)) {
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) vecval[0] * one[0];
		    d2 = (double) -vecval[1] * one[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) vecval[0] * one[1];
		    d2 = (double) vecval[1] * one[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) matval[0] * vecval[0];
		    d2 = (double) -matval[1] * vecval[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) matval[0] * vecval[1];
		    d2 = (double) matval[1] * vecval[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_rowsum[0];
		  tail_a = tail_rowsum[0];
		  head_b = head_rowtmp[0];
		  tail_b = tail_rowtmp[0];
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
		  head_rowsum[0] = head_t;
		  tail_rowsum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_rowsum[1];
		  tail_a = tail_rowsum[1];
		  head_b = head_rowtmp[1];
		  tail_b = tail_rowtmp[1];
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
		  head_rowsum[1] = head_t;
		  tail_rowsum[1] = tail_t;
		}
		x_index2 -= incx;
		tp_index -= inctp;
	      }
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
		  head_a0 = head_rowsum[0];
		  tail_a0 = tail_rowsum[0];
		  head_a1 = head_rowsum[1];
		  tail_a1 = tail_rowsum[1];
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
		  head_result[0] = head_t1;
		  tail_result[0] = tail_t1;
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
		  head_result[1] = head_t1;
		  tail_result[1] = tail_t1;
		}

	      }
	      x_i[x_index] = head_result[0];
	      x_i[x_index + 1] = head_result[1];
	      x_index -= incx;
	    }
	  } else {
	    tp_start = 0;
	    x_index = x_start + (n - 1) * incx;
	    for (matrix_row = n - 1; matrix_row >= 0; matrix_row--) {
	      tp_index = matrix_row * inctp;
	      x_index2 = x_start;
	      head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		tail_rowsum[1] = 0.0;
	      head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		tail_rowtmp[1] = 0.0;
	      head_result[0] = head_result[1] = tail_result[0] =
		tail_result[1] = 0.0;
	      stride = n;
	      for (step = 0; step <= matrix_row; step++) {
		vecval[0] = x_i[x_index2];
		vecval[1] = x_i[x_index2 + 1];
		if ((diag == blas_unit_diag) && (step == matrix_row)) {
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) vecval[0] * one[0];
		    d2 = (double) -vecval[1] * one[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) vecval[0] * one[1];
		    d2 = (double) vecval[1] * one[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		} else {
		  matval[0] = tp_i[tp_index];
		  matval[1] = tp_i[tp_index + 1];
		  {
		    double head_e1, tail_e1;
		    double d1;
		    double d2;
		    /* Real part */
		    d1 = (double) matval[0] * vecval[0];
		    d2 = (double) -matval[1] * vecval[1];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[0] = head_e1;
		    tail_rowtmp[0] = tail_e1;
		    /* imaginary part */
		    d1 = (double) matval[0] * vecval[1];
		    d2 = (double) matval[1] * vecval[0];
		    {
		      /* Compute double-double = double + double. */
		      double e, t1, t2;

		      /* Knuth trick. */
		      t1 = d1 + d2;
		      e = t1 - d1;
		      t2 = ((d2 - e) + (d1 - (t1 - e)));

		      /* The result is t1 + t2, after normalization. */
		      head_e1 = t1 + t2;
		      tail_e1 = t2 - (head_e1 - t1);
		    }
		    head_rowtmp[1] = head_e1;
		    tail_rowtmp[1] = tail_e1;
		  }
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_rowsum[0];
		  tail_a = tail_rowsum[0];
		  head_b = head_rowtmp[0];
		  tail_b = tail_rowtmp[0];
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
		  head_rowsum[0] = head_t;
		  tail_rowsum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_rowsum[1];
		  tail_a = tail_rowsum[1];
		  head_b = head_rowtmp[1];
		  tail_b = tail_rowtmp[1];
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
		  head_rowsum[1] = head_t;
		  tail_rowsum[1] = tail_t;
		}
		stride--;
		tp_index += stride * inctp;
		x_index2 += incx;
	      }
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
		  head_a0 = head_rowsum[0];
		  tail_a0 = tail_rowsum[0];
		  head_a1 = head_rowsum[1];
		  tail_a1 = tail_rowsum[1];
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
		  head_result[0] = head_t1;
		  tail_result[0] = tail_t1;
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
		  head_result[1] = head_t1;
		  tail_result[1] = tail_t1;
		}

	      }
	      x_i[x_index] = head_result[0];
	      x_i[x_index + 1] = head_result[1];
	      x_index -= incx;
	    }
	  }
	}


      }
      break;
    }
  }

}
