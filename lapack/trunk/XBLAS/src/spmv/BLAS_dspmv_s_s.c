#include "blas_extended.h"
#include "blas_extended_private.h"

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * ap * x + beta * y, where ap is a symmetric
 * packed matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of ap; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether ap is upper or lower
 *
 * n            (input) int
 *              Dimension of ap and the length of vector x
 *
 * alpha        (input) double
 *              
 * ap           (input) float*
 *
 * x            (input) float*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) double
 *
 * y            (input/output) double*
 *
 * incy         (input) int
 *              The stride for vector y.
 *
 */
void BLAS_dspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *ap,
		    const float *x, int incx, double beta,
		    double *y, int incy)
{
  static const char routine_name[] = "BLAS_dspmv_s_s";

  {
    int matrix_row, step, ap_index, ap_start, x_index, x_start;
    int y_start, y_index, incap;
    double alpha_i = alpha;
    double beta_i = beta;

    const float *ap_i = ap;
    const float *x_i = x;
    double *y_i = y;
    double rowsum;
    double rowtmp;
    float matval;
    float vecval;
    double resval;
    double tmp1;
    double tmp2;


    incap = 1;




    if (incx < 0)
      x_start = (-n + 1) * incx;
    else
      x_start = 0;
    if (incy < 0)
      y_start = (-n + 1) * incy;
    else
      y_start = 0;

    if (n < 1) {
      return;
    }
    if (alpha_i == 0.0 && beta_i == 1.0) {
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
      BLAS_error(routine_name, -7, incx, NULL);
    }
    if (incy == 0) {
      BLAS_error(routine_name, -10, incy, NULL);
    }



    if (alpha_i == 0.0) {
      {
	y_index = y_start;
	for (matrix_row = 0; matrix_row < n; matrix_row++) {
	  resval = y_i[y_index];

	  tmp2 = beta_i * resval;

	  y_i[y_index] = tmp2;

	  y_index += incy;
	}
      }
    } else {
      if (uplo == blas_lower)
	order = (order == blas_rowmajor) ? blas_colmajor : blas_rowmajor;
      if (order == blas_rowmajor) {
	if (alpha_i == 1.0) {
	  if (beta_i == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		tmp1 = rowsum;
		y_i[y_index] = tmp1;

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  } else {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		resval = y_i[y_index];
		tmp1 = rowsum;
		tmp2 = beta_i * resval;
		tmp2 = tmp1 + tmp2;
		y_i[y_index] = tmp2;

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  }
	} else {
	  if (beta_i == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		tmp1 = rowsum * alpha_i;
		y_i[y_index] = tmp1;

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  } else {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		resval = y_i[y_index];
		tmp1 = rowsum * alpha_i;
		tmp2 = beta_i * resval;
		tmp2 = tmp1 + tmp2;
		y_i[y_index] = tmp2;

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  }
	}
      } else {
	if (alpha_i == 1.0) {
	  if (beta_i == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		tmp1 = rowsum;
		y_i[y_index] = tmp1;

		y_index += incy;
		ap_start += (matrix_row + 1) * incap;
	      }
	    }
	  } else {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		resval = y_i[y_index];
		tmp1 = rowsum;
		tmp2 = beta_i * resval;
		tmp2 = tmp1 + tmp2;
		y_i[y_index] = tmp2;

		y_index += incy;
		ap_start += (matrix_row + 1) * incap;
	      }
	    }
	  }
	} else {
	  if (beta_i == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		tmp1 = rowsum * alpha_i;
		y_i[y_index] = tmp1;

		y_index += incy;
		ap_start += (matrix_row + 1) * incap;
	      }
	    }
	  } else {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum = 0.0;
		rowtmp = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval = ap_i[ap_index];
		  vecval = x_i[x_index];
		  rowtmp = (double) matval *vecval;
		  rowsum = rowsum + rowtmp;
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		resval = y_i[y_index];
		tmp1 = rowsum * alpha_i;
		tmp2 = beta_i * resval;
		tmp2 = tmp1 + tmp2;
		y_i[y_index] = tmp2;

		y_index += incy;
		ap_start += (matrix_row + 1) * incap;
	      }
	    }
	  }
	}
      }				/* if order == ... */
    }				/* alpha != 0 */


  }
}
