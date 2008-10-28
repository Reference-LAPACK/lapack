#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_ctrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const float *T, int ldt,
		  void *x, int incx)

/*
 * Purpose
 * =======
 * 
 * This routine solve :
 * 
 *     x <- alpha * inverse(T) * x
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
 * x      (input) const void*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 */
{
  char *routine_name = "BLAS_ctrsv_s";

  int i, j;			/* used to idx matrix */
  int ix, jx;			/* used to idx vector x */
  int start_x;			/* used as the starting idx to vector x */
  const float *T_i = T;		/* internal matrix T */
  float *x_i = (float *) x;	/* internal x */
  float *alpha_i = (float *) alpha;	/* internal alpha */
  float T_element;		/* temporary variable for an element of matrix A */
  int incT = 1;			/* internal ldt */

  if ((order != blas_rowmajor && order != blas_colmajor) ||
      (uplo != blas_upper && uplo != blas_lower) ||
      (trans != blas_trans && trans !=
       blas_no_trans && trans != blas_conj_trans) ||
      (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
      (ldt < n) || (incx == 0)) {
    BLAS_error(routine_name, 0, 0, NULL);
  }

  if (n <= 0)
    return;


  incx *= 2;
  /* configuring the vector starting idx */
  if (incx <= 0) {
    start_x = -(n - 1) * incx;
  } else {
    start_x = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    ix = start_x;
    for (i = 0; i < n; i++) {
      x_i[ix] = 0.0;
      x_i[ix + 1] = 0.0;
      ix += incx;
    }
    return;
  }

  {
    float temp1[2];		/* temporary variable for calculations */
    float temp2[2];		/* temporary variable for calculations */
    float temp3[2];		/* temporary variable for calculations */

    if ((order == blas_rowmajor &&
	 trans == blas_no_trans && uplo == blas_upper) ||
	(order == blas_colmajor &&
	 trans != blas_no_trans && uplo == blas_lower)) {

      jx = start_x + (n - 1) * incx;
      for (j = n - 1; j >= 0; j--) {

	/* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	   i=j+1 to n-1           */
	temp3[0] = x_i[jx];
	temp3[1] = x_i[1 + jx];
	{
	  temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	  temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	}


	ix = start_x + (n - 1) * incx;
	for (i = n - 1; i >= j + 1; i--) {
	  T_element = T_i[i * incT + j * ldt * incT];

	  temp3[0] = x_i[ix];
	  temp3[1] = x_i[1 + ix];
	  {
	    temp2[0] = temp3[0] * T_element;
	    temp2[1] = temp3[1] * T_element;
	  }
	  temp1[0] = temp1[0] - temp2[0];
	  temp1[1] = temp1[1] - temp2[1];
	  ix -= incx;
	}			/* for j<n */

	/* if the diagonal entry is not equal to one, then divide Xj by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = T_i[j * incT + j * ldt * incT];


	  temp1[0] = temp1[0] / T_element;
	  temp1[1] = temp1[1] / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[jx] = temp1[0];
	x_i[jx + 1] = temp1[1];

	jx -= incx;
      }				/* for j>=0 */
    } else if ((order == blas_rowmajor &&
		trans == blas_no_trans && uplo == blas_lower) ||
	       (order == blas_colmajor &&
		trans != blas_no_trans && uplo == blas_upper)) {

      jx = start_x;
      for (j = 0; j < n; j++) {

	/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	   i=j+1 to n-1           */
	temp3[0] = x_i[jx];
	temp3[1] = x_i[1 + jx];
	/* multiply by alpha */
	{
	  temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	  temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	}


	ix = start_x;
	for (i = 0; i < j; i++) {
	  T_element = T_i[i * incT + j * ldt * incT];

	  temp3[0] = x_i[ix];
	  temp3[1] = x_i[1 + ix];
	  {
	    temp2[0] = temp3[0] * T_element;
	    temp2[1] = temp3[1] * T_element;
	  }
	  temp1[0] = temp1[0] - temp2[0];
	  temp1[1] = temp1[1] - temp2[1];
	  ix += incx;
	}			/* for i<j */

	/* if the diagonal entry is not equal to one, then divide Xj by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = T_i[j * incT + j * ldt * incT];


	  temp1[0] = temp1[0] / T_element;
	  temp1[1] = temp1[1] / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[jx] = temp1[0];
	x_i[jx + 1] = temp1[1];
	jx += incx;
      }				/* for j<n */
    } else if ((order == blas_rowmajor &&
		trans != blas_no_trans && uplo == blas_lower) ||
	       (order == blas_colmajor &&
		trans == blas_no_trans && uplo == blas_upper)) {

      jx = start_x + (n - 1) * incx;
      for (j = n - 1; j >= 0; j--) {

	/* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	   i=j+1 to n-1           */
	temp3[0] = x_i[jx];
	temp3[1] = x_i[1 + jx];
	{
	  temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	  temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	}


	ix = start_x + (n - 1) * incx;
	for (i = n - 1; i >= j + 1; i--) {
	  T_element = T_i[j * incT + i * ldt * incT];

	  temp3[0] = x_i[ix];
	  temp3[1] = x_i[1 + ix];
	  {
	    temp2[0] = temp3[0] * T_element;
	    temp2[1] = temp3[1] * T_element;
	  }
	  temp1[0] = temp1[0] - temp2[0];
	  temp1[1] = temp1[1] - temp2[1];
	  ix -= incx;
	}			/* for j<n */

	/* if the diagonal entry is not equal to one, then divide Xj by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = T_i[j * incT + j * ldt * incT];


	  temp1[0] = temp1[0] / T_element;
	  temp1[1] = temp1[1] / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[jx] = temp1[0];
	x_i[jx + 1] = temp1[1];

	jx -= incx;
      }				/* for j>=0 */
    } else if ((order == blas_rowmajor &&
		trans != blas_no_trans && uplo == blas_upper) ||
	       (order == blas_colmajor &&
		trans == blas_no_trans && uplo == blas_lower)) {

      jx = start_x;
      for (j = 0; j < n; j++) {

	/* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	   i=j+1 to n-1           */
	temp3[0] = x_i[jx];
	temp3[1] = x_i[1 + jx];
	/* multiply by alpha */
	{
	  temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	  temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	}


	ix = start_x;
	for (i = 0; i < j; i++) {
	  T_element = T_i[j * incT + i * ldt * incT];

	  temp3[0] = x_i[ix];
	  temp3[1] = x_i[1 + ix];
	  {
	    temp2[0] = temp3[0] * T_element;
	    temp2[1] = temp3[1] * T_element;
	  }
	  temp1[0] = temp1[0] - temp2[0];
	  temp1[1] = temp1[1] - temp2[1];
	  ix += incx;
	}			/* for i<j */

	/* if the diagonal entry is not equal to one, then divide Xj by 
	   the entry */
	if (diag == blas_non_unit_diag) {
	  T_element = T_i[j * incT + j * ldt * incT];


	  temp1[0] = temp1[0] / T_element;
	  temp1[1] = temp1[1] / T_element;

	}
	/* if (diag == blas_non_unit_diag) */
	x_i[jx] = temp1[0];
	x_i[jx + 1] = temp1[1];
	jx += incx;
      }				/* for j<n */
    }
  }
}
