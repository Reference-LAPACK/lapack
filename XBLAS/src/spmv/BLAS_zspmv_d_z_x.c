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
 * ap           (input) double*
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
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 *
 */
void BLAS_zspmv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec)
{
  static const char routine_name[] = "BLAS_zspmv_d_z_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_indigenous:
  case blas_prec_double:{
      {
	int matrix_row, step, ap_index, ap_start, x_index, x_start;
	int y_start, y_index, incap;
	double *alpha_i = (double *) alpha;
	double *beta_i = (double *) beta;

	const double *ap_i = ap;
	const double *x_i = (double *) x;
	double *y_i = (double *) y;
	double rowsum[2];
	double rowtmp[2];
	double matval;
	double vecval[2];
	double resval[2];
	double tmp1[2];
	double tmp2[2];


	incap = 1;

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
		  (double) beta_i[0] * resval[0] -
		  (double) beta_i[1] * resval[1];
		tmp2[1] =
		  (double) beta_i[0] * resval[1] +
		  (double) beta_i[1] * resval[0];
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
		      }
		      rowsum[0] = rowsum[0] + rowtmp[0];
		      rowsum[1] = rowsum[1] + rowtmp[1];
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			rowtmp[0] = vecval[0] * matval;
			rowtmp[1] = vecval[1] * matval;
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
	  }			/* if order == ... */
	}			/* alpha != 0 */


      }
      break;
    }

  case blas_prec_extra:{
      {
	int matrix_row, step, ap_index, ap_start, x_index, x_start;
	int y_start, y_index, incap;
	double *alpha_i = (double *) alpha;
	double *beta_i = (double *) beta;

	const double *ap_i = ap;
	const double *x_i = (double *) x;
	double *y_i = (double *) y;
	double head_rowsum[2], tail_rowsum[2];
	double head_rowtmp[2], tail_rowtmp[2];
	double matval;
	double vecval[2];
	double resval[2];
	double head_tmp1[2], tail_tmp1[2];
	double head_tmp2[2], tail_tmp2[2];
	FPU_FIX_DECL;

	incap = 1;

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

	FPU_FIX_START;

	if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	  {
	    y_index = y_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      resval[0] = y_i[y_index];
	      resval[1] = y_i[y_index + 1];

	      {
		/* Compute complex-extra = complex-double * complex-double. */
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		/* Real part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = beta_i[0] * split;
		  a1 = con - beta_i[0];
		  a1 = con - a1;
		  a2 = beta_i[0] - a1;
		  con = resval[0] * split;
		  b1 = con - resval[0];
		  b1 = con - b1;
		  b2 = resval[0] - b1;

		  head_t1 = beta_i[0] * resval[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = beta_i[1] * split;
		  a1 = con - beta_i[1];
		  a1 = con - a1;
		  a2 = beta_i[1] - a1;
		  con = resval[1] * split;
		  b1 = con - resval[1];
		  b1 = con - b1;
		  b2 = resval[1] - b1;

		  head_t2 = beta_i[1] * resval[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_tmp2[0] = head_t1;
		tail_tmp2[0] = tail_t1;
		/* Imaginary part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = beta_i[1] * split;
		  a1 = con - beta_i[1];
		  a1 = con - a1;
		  a2 = beta_i[1] - a1;
		  con = resval[0] * split;
		  b1 = con - resval[0];
		  b1 = con - b1;
		  b2 = resval[0] - b1;

		  head_t1 = beta_i[1] * resval[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = beta_i[0] * split;
		  a1 = con - beta_i[0];
		  a1 = con - a1;
		  a2 = beta_i[0] - a1;
		  con = resval[1] * split;
		  b1 = con - resval[1];
		  b1 = con - b1;
		  b2 = resval[1] - b1;

		  head_t2 = beta_i[0] * resval[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		head_tmp2[1] = head_t1;
		tail_tmp2[1] = tail_t1;
	      }

	      y_i[y_index] = head_tmp2[0];
	      y_i[y_index + 1] = head_tmp2[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    head_tmp1[0] = head_rowsum[0];
		    tail_tmp1[0] = tail_rowsum[0];
		    head_tmp1[1] = head_rowsum[1];
		    tail_tmp1[1] = tail_rowsum[1];
		    y_i[y_index] = head_tmp1[0];
		    y_i[y_index + 1] = head_tmp1[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval[0] = y_i[y_index];
		    resval[1] = y_i[y_index + 1];
		    head_tmp1[0] = head_rowsum[0];
		    tail_tmp1[0] = tail_rowsum[0];
		    head_tmp1[1] = head_rowsum[1];
		    tail_tmp1[1] = tail_rowsum[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[0] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[1] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[0] = head_t1;
		      tail_tmp2[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[1] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[0] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[1] = head_t1;
		      tail_tmp2[1] = tail_t1;
		    }
		    {
		      double head_t, tail_t;
		      double head_a, tail_a;
		      double head_b, tail_b;
		      /* Real part */
		      head_a = head_tmp1[0];
		      tail_a = tail_tmp1[0];
		      head_b = head_tmp2[0];
		      tail_b = tail_tmp2[0];
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
		      head_tmp2[0] = head_t;
		      tail_tmp2[0] = tail_t;
		      /* Imaginary part */
		      head_a = head_tmp1[1];
		      tail_a = tail_tmp1[1];
		      head_b = head_tmp2[1];
		      tail_b = tail_tmp2[1];
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
		      head_tmp2[1] = head_t;
		      tail_tmp2[1] = tail_t;
		    }
		    y_i[y_index] = head_tmp2[0];
		    y_i[y_index + 1] = head_tmp2[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
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
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a0 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a1 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[1];
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
		      head_tmp1[0] = head_t1;
		      tail_tmp1[0] = tail_t1;
		      /* imaginary part */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a1 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a0 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[1];
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
		      head_tmp1[1] = head_t1;
		      tail_tmp1[1] = tail_t1;
		    }

		    y_i[y_index] = head_tmp1[0];
		    y_i[y_index + 1] = head_tmp1[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval[0] = y_i[y_index];
		    resval[1] = y_i[y_index + 1];
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
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a0 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a1 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[1];
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
		      head_tmp1[0] = head_t1;
		      tail_tmp1[0] = tail_t1;
		      /* imaginary part */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a1 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a0 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[1];
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
		      head_tmp1[1] = head_t1;
		      tail_tmp1[1] = tail_t1;
		    }

		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[0] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[1] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[0] = head_t1;
		      tail_tmp2[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[1] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[0] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[1] = head_t1;
		      tail_tmp2[1] = tail_t1;
		    }
		    {
		      double head_t, tail_t;
		      double head_a, tail_a;
		      double head_b, tail_b;
		      /* Real part */
		      head_a = head_tmp1[0];
		      tail_a = tail_tmp1[0];
		      head_b = head_tmp2[0];
		      tail_b = tail_tmp2[0];
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
		      head_tmp2[0] = head_t;
		      tail_tmp2[0] = tail_t;
		      /* Imaginary part */
		      head_a = head_tmp1[1];
		      tail_a = tail_tmp1[1];
		      head_b = head_tmp2[1];
		      tail_b = tail_tmp2[1];
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
		      head_tmp2[1] = head_t;
		      tail_tmp2[1] = tail_t;
		    }
		    y_i[y_index] = head_tmp2[0];
		    y_i[y_index + 1] = head_tmp2[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    head_tmp1[0] = head_rowsum[0];
		    tail_tmp1[0] = tail_rowsum[0];
		    head_tmp1[1] = head_rowsum[1];
		    tail_tmp1[1] = tail_rowsum[1];
		    y_i[y_index] = head_tmp1[0];
		    y_i[y_index + 1] = head_tmp1[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval[0] = y_i[y_index];
		    resval[1] = y_i[y_index + 1];
		    head_tmp1[0] = head_rowsum[0];
		    tail_tmp1[0] = tail_rowsum[0];
		    head_tmp1[1] = head_rowsum[1];
		    tail_tmp1[1] = tail_rowsum[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[0] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[1] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[0] = head_t1;
		      tail_tmp2[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[1] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[0] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[1] = head_t1;
		      tail_tmp2[1] = tail_t1;
		    }
		    {
		      double head_t, tail_t;
		      double head_a, tail_a;
		      double head_b, tail_b;
		      /* Real part */
		      head_a = head_tmp1[0];
		      tail_a = tail_tmp1[0];
		      head_b = head_tmp2[0];
		      tail_b = tail_tmp2[0];
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
		      head_tmp2[0] = head_t;
		      tail_tmp2[0] = tail_t;
		      /* Imaginary part */
		      head_a = head_tmp1[1];
		      tail_a = tail_tmp1[1];
		      head_b = head_tmp2[1];
		      tail_b = tail_tmp2[1];
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
		      head_tmp2[1] = head_t;
		      tail_tmp2[1] = tail_t;
		    }
		    y_i[y_index] = head_tmp2[0];
		    y_i[y_index + 1] = head_tmp2[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
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
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a0 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a1 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[1];
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
		      head_tmp1[0] = head_t1;
		      tail_tmp1[0] = tail_t1;
		      /* imaginary part */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a1 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a0 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[1];
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
		      head_tmp1[1] = head_t1;
		      tail_tmp1[1] = tail_t1;
		    }

		    y_i[y_index] = head_tmp1[0];
		    y_i[y_index + 1] = head_tmp1[1];

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
		    head_rowsum[0] = head_rowsum[1] = tail_rowsum[0] =
		      tail_rowsum[1] = 0.0;
		    head_rowtmp[0] = head_rowtmp[1] = tail_rowtmp[0] =
		      tail_rowtmp[1] = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval[0] = x_i[x_index];
		      vecval[1] = x_i[x_index + 1];
		      {
			/* Compute complex-extra = complex-double * real. */
			double head_t, tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[0] * split;
			  b1 = con - vecval[0];
			  b1 = con - b1;
			  b2 = vecval[0] - b1;

			  head_t = matval * vecval[0];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[0] = head_t;
			tail_rowtmp[0] = tail_t;
			{
			  /* Compute double_double = double * double. */
			  double a1, a2, b1, b2, con;

			  con = matval * split;
			  a1 = con - matval;
			  a1 = con - a1;
			  a2 = matval - a1;
			  con = vecval[1] * split;
			  b1 = con - vecval[1];
			  b1 = con - b1;
			  b2 = vecval[1] - b1;

			  head_t = matval * vecval[1];
			  tail_t =
			    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) +
			    a2 * b2;
			}
			head_rowtmp[1] = head_t;
			tail_rowtmp[1] = tail_t;
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
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval[0] = y_i[y_index];
		    resval[1] = y_i[y_index + 1];
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
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a0 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a1 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[1];
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
		      head_tmp1[0] = head_t1;
		      tail_tmp1[0] = tail_t1;
		      /* imaginary part */
		      {
			/* Compute double-double = double-double * double. */
			double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

			con = head_a1 * split;
			a11 = con - head_a1;
			a11 = con - a11;
			a21 = head_a1 - a11;
			con = alpha_i[0] * split;
			b1 = con - alpha_i[0];
			b1 = con - b1;
			b2 = alpha_i[0] - b1;

			c11 = head_a1 * alpha_i[0];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a1 * alpha_i[0];
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
			con = alpha_i[1] * split;
			b1 = con - alpha_i[1];
			b1 = con - b1;
			b2 = alpha_i[1] - b1;

			c11 = head_a0 * alpha_i[1];
			c21 =
			  (((a11 * b1 - c11) + a11 * b2) + a21 * b1) +
			  a21 * b2;

			c2 = tail_a0 * alpha_i[1];
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
		      head_tmp1[1] = head_t1;
		      tail_tmp1[1] = tail_t1;
		    }

		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[0] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[1] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[0] = head_t1;
		      tail_tmp2[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[1] * split;
			a1 = con - beta_i[1];
			a1 = con - a1;
			a2 = beta_i[1] - a1;
			con = resval[0] * split;
			b1 = con - resval[0];
			b1 = con - b1;
			b2 = resval[0] - b1;

			head_t1 = beta_i[1] * resval[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = beta_i[0] * split;
			a1 = con - beta_i[0];
			a1 = con - a1;
			a2 = beta_i[0] - a1;
			con = resval[1] * split;
			b1 = con - resval[1];
			b1 = con - b1;
			b2 = resval[1] - b1;

			head_t2 = beta_i[0] * resval[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_tmp2[1] = head_t1;
		      tail_tmp2[1] = tail_t1;
		    }
		    {
		      double head_t, tail_t;
		      double head_a, tail_a;
		      double head_b, tail_b;
		      /* Real part */
		      head_a = head_tmp1[0];
		      tail_a = tail_tmp1[0];
		      head_b = head_tmp2[0];
		      tail_b = tail_tmp2[0];
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
		      head_tmp2[0] = head_t;
		      tail_tmp2[0] = tail_t;
		      /* Imaginary part */
		      head_a = head_tmp1[1];
		      tail_a = tail_tmp1[1];
		      head_b = head_tmp2[1];
		      tail_b = tail_tmp2[1];
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
		      head_tmp2[1] = head_t;
		      tail_tmp2[1] = tail_t;
		    }
		    y_i[y_index] = head_tmp2[0];
		    y_i[y_index + 1] = head_tmp2[1];

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      }
	    }
	  }			/* if order == ... */
	}			/* alpha != 0 */

	FPU_FIX_STOP;
      }
      break;
    }

  }
}
