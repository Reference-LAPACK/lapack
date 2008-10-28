#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_dtpmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const float *tp, double *x, int incx)

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
 * alpha        (input) double
 *
 * tp           (input) float*
 *
 * x            (input) double*
 *
 * incx         (input) int
 *              The stride for vector x.
 *
 */
{
  static const char routine_name[] = "BLAS_dtpmv_s";

  {
    int matrix_row, step, tp_index, tp_start, x_index, x_start;
    int inctp, x_index2, stride, col_index, inctp2;

    double alpha_i = alpha;

    const float *tp_i = tp;
    double *x_i = x;
    double rowsum;
    double rowtmp;
    double result;
    float matval;
    double vecval;
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
}
