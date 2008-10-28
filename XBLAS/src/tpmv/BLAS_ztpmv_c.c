#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ztpmv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp, void *x, int incx)

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
 */
{
  static const char routine_name[] = "BLAS_ztpmv_c";

  {
    int matrix_row, step, tp_index, tp_start, x_index, x_start;
    int inctp, x_index2, stride, col_index, inctp2;

    double *alpha_i = (double *) alpha;

    const float *tp_i = (float *) tp;
    double *x_i = (double *) x;
    double rowsum[2];
    double rowtmp[2];
    double result[2];
    float matval[2];
    double vecval[2];
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
		  (double) vecval[0] * one[0] - (double) vecval[1] * one[1];
		rowtmp[1] =
		  (double) vecval[0] * one[1] + (double) vecval[1] * one[0];
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
		  (double) vecval[0] * one[0] - (double) vecval[1] * one[1];
		rowtmp[1] =
		  (double) vecval[0] * one[1] + (double) vecval[1] * one[0];
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
		  (double) vecval[0] * one[0] - (double) vecval[1] * one[1];
		rowtmp[1] =
		  (double) vecval[0] * one[1] + (double) vecval[1] * one[0];
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
		  (double) vecval[0] * one[0] - (double) vecval[1] * one[1];
		rowtmp[1] =
		  (double) vecval[0] * one[1] + (double) vecval[1] * one[0];
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
}
