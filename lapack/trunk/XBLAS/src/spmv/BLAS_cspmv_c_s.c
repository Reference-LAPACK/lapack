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
 * x            (input) float*
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
void BLAS_cspmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const float *x, int incx, const void *beta,
		    void *y, int incy)
{
  static const char routine_name[] = "BLAS_cspmv_c_s";

  {
    int matrix_row, step, ap_index, ap_start, x_index, x_start;
    int y_start, y_index, incap;
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;

    const float *ap_i = (float *) ap;
    const float *x_i = x;
    float *y_i = (float *) y;
    float rowsum[2];
    float rowtmp[2];
    float matval[2];
    float vecval;
    float resval[2];
    float tmp1[2];
    float tmp2[2];


    incap = 1;
    incap *= 2;

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
	    tmp2[0] = beta_i[0] * resval[0] - beta_i[1] * resval[1];
	    tmp2[1] = beta_i[0] * resval[1] + beta_i[1] * resval[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
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
		  tmp2[0] = beta_i[0] * resval[0] - beta_i[1] * resval[1];
		  tmp2[1] = beta_i[0] * resval[1] + beta_i[1] * resval[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		{
		  tmp1[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		  tmp1[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (n - step - 1) * incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		{
		  tmp1[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		  tmp1[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
		}

		{
		  tmp2[0] = beta_i[0] * resval[0] - beta_i[1] * resval[1];
		  tmp2[1] = beta_i[0] * resval[1] + beta_i[1] * resval[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
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
		  tmp2[0] = beta_i[0] * resval[0] - beta_i[1] * resval[1];
		  tmp2[1] = beta_i[0] * resval[1] + beta_i[1] * resval[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		{
		  tmp1[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		  tmp1[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
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
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += incap;
		  x_index += incx;
		}
		for (step = matrix_row; step < n; step++) {
		  matval[0] = ap_i[ap_index];
		  matval[1] = ap_i[ap_index + 1];
		  vecval = x_i[x_index];
		  {
		    rowtmp[0] = matval[0] * vecval;
		    rowtmp[1] = matval[1] * vecval;
		  }
		  rowsum[0] = rowsum[0] + rowtmp[0];
		  rowsum[1] = rowsum[1] + rowtmp[1];
		  ap_index += (step + 1) * incap;
		  x_index += incx;
		}
		resval[0] = y_i[y_index];
		resval[1] = y_i[y_index + 1];
		{
		  tmp1[0] = rowsum[0] * alpha_i[0] - rowsum[1] * alpha_i[1];
		  tmp1[1] = rowsum[0] * alpha_i[1] + rowsum[1] * alpha_i[0];
		}

		{
		  tmp2[0] = beta_i[0] * resval[0] - beta_i[1] * resval[1];
		  tmp2[1] = beta_i[0] * resval[1] + beta_i[1] * resval[0];
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
