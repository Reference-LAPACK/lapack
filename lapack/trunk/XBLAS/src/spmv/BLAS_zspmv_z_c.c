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
 * alpha        (input) const void*
 *              
 * ap           (input) void*
 *
 * x            (input) void*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) const void*
 *
 * y            (input/output) void*
 *
 * incy         (input) int
 *              The stride for vector y.
 *
 */
void BLAS_zspmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy)
{
  static const char routine_name[] = "BLAS_zspmv_z_c";

  {
    int matrix_row, step, ap_index, ap_start, x_index, x_start;
    int y_start, y_index, incap;
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;

    const double *ap_i = (double *) ap;
    const float *x_i = (float *) x;
    double *y_i = (double *) y;
    double rowsum[2];
    double rowtmp[2];
    double matval[2];
    float vecval[2];
    double resval[2];
    double tmp1[2];
    double tmp2[2];


    incap = 1;
    incap *= 2;
    incx *= 2;
    incy *= 2;

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
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	&& (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
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



    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      {
	y_index = y_start;
	for (matrix_row = 0; matrix_row < n; matrix_row++) {
	  resval[0] = y_i[y_index];
	  resval[1] = y_i[y_index + 1];

	  {
	    tmp2[0] =
	      (double) beta_i[0] * resval[0] - (double) beta_i[1] * resval[1];
	    tmp2[1] =
	      (double) beta_i[0] * resval[1] + (double) beta_i[1] * resval[0];
	  }

	  y_i[y_index] = tmp2[0];
	  y_i[y_index + 1] = tmp2[1];

	  y_index += incy;
	}
      }
    } else {
      if (uplo == blas_lower)
	order = (order == blas_rowmajor) ? blas_colmajor : blas_rowmajor;
      if (order == blas_rowmajor) {
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		tmp1[0] = rowsum[0];
		tmp1[1] = rowsum[1];
		y_i[y_index] = tmp1[0];
		y_i[y_index + 1] = tmp1[1];

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
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		tmp1[0] = rowsum[0];
		tmp1[1] = rowsum[1];
		{
		  tmp2[0] =
		    (double) beta_i[0] * resval[0] -
		    (double) beta_i[1] * resval[1];
		  tmp2[1] =
		    (double) beta_i[0] * resval[1] +
		    (double) beta_i[1] * resval[0];
		}
		tmp2[0] = tmp1[0] + tmp2[0];
		tmp2[1] = tmp1[1] + tmp2[1];
		y_i[y_index] = tmp2[0];
		y_i[y_index + 1] = tmp2[1];

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  }
	} else {
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		{
		  tmp1[0] =
		    (double) rowsum[0] * alpha_i[0] -
		    (double) rowsum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) rowsum[0] * alpha_i[1] +
		    (double) rowsum[1] * alpha_i[0];
		}
		y_i[y_index] = tmp1[0];
		y_i[y_index + 1] = tmp1[1];

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
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		{
		  tmp1[0] =
		    (double) rowsum[0] * alpha_i[0] -
		    (double) rowsum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) rowsum[0] * alpha_i[1] +
		    (double) rowsum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) beta_i[0] * resval[0] -
		    (double) beta_i[1] * resval[1];
		  tmp2[1] =
		    (double) beta_i[0] * resval[1] +
		    (double) beta_i[1] * resval[0];
		}
		tmp2[0] = tmp1[0] + tmp2[0];
		tmp2[1] = tmp1[1] + tmp2[1];
		y_i[y_index] = tmp2[0];
		y_i[y_index + 1] = tmp2[1];

		y_index += incy;
		ap_start += incap;
	      }
	    }
	  }
	}
      } else {
	if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		tmp1[0] = rowsum[0];
		tmp1[1] = rowsum[1];
		y_i[y_index] = tmp1[0];
		y_i[y_index + 1] = tmp1[1];

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
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		tmp1[0] = rowsum[0];
		tmp1[1] = rowsum[1];
		{
		  tmp2[0] =
		    (double) beta_i[0] * resval[0] -
		    (double) beta_i[1] * resval[1];
		  tmp2[1] =
		    (double) beta_i[0] * resval[1] +
		    (double) beta_i[1] * resval[0];
		}
		tmp2[0] = tmp1[0] + tmp2[0];
		tmp2[1] = tmp1[1] + tmp2[1];
		y_i[y_index] = tmp2[0];
		y_i[y_index + 1] = tmp2[1];

		y_index += incy;
		ap_start += (matrix_row + 1) * incap;
	      }
	    }
	  }
	} else {
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    {
	      y_index = y_start;
	      ap_start = 0;
	      for (matrix_row = 0; matrix_row < n; matrix_row++) {
		x_index = x_start;
		ap_index = ap_start;
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		{
		  tmp1[0] =
		    (double) rowsum[0] * alpha_i[0] -
		    (double) rowsum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) rowsum[0] * alpha_i[1] +
		    (double) rowsum[1] * alpha_i[0];
		}
		y_i[y_index] = tmp1[0];
		y_i[y_index + 1] = tmp1[1];

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
		rowsum[0] = rowsum[1] = 0.0;
		rowtmp[0] = rowtmp[1] = 0.0;
		for (step = 0; step < matrix_row; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval[0] = x_i[x_index];
		  vecval[1] = x_i[x_index + 1];
		  {
		    rowtmp[0] =
		      (double) matval[0] * vecval[0] -
		      (double) matval[1] * vecval[1];
		    rowtmp[1] =
		      (double) matval[0] * vecval[1] +
		      (double) matval[1] * vecval[0];
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		{
		  tmp1[0] =
		    (double) rowsum[0] * alpha_i[0] -
		    (double) rowsum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) rowsum[0] * alpha_i[1] +
		    (double) rowsum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) beta_i[0] * resval[0] -
		    (double) beta_i[1] * resval[1];
		  tmp2[1] =
		    (double) beta_i[0] * resval[1] +
		    (double) beta_i[1] * resval[0];
		}
		tmp2[0] = tmp1[0] + tmp2[0];
		tmp2[1] = tmp1[1] + tmp2[1];
		y_i[y_index] = tmp2[0];
		y_i[y_index + 1] = tmp2[1];

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
