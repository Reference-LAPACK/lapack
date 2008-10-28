#include "blas_extended.h"
#include "blas_extended_test.h"

void strsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int length, float *T, int lda,
		  const float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i;
  float tmp;

  if ((order == blas_rowmajor && uplo == blas_upper && trans == blas_no_trans)
      || (order == blas_colmajor && uplo == blas_lower
	  && trans != blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row * lda) + (row + i + 1)] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_lower
	 && trans == blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_upper
					&& trans != blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_lower
	 && trans != blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_upper
					&& trans == blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_upper
	 && trans != blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_lower
					&& trans == blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  }
}

void dtrsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int length, double *T, int lda,
		  const double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i;
  double tmp;

  if ((order == blas_rowmajor && uplo == blas_upper && trans == blas_no_trans)
      || (order == blas_colmajor && uplo == blas_lower
	  && trans != blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row * lda) + (row + i + 1)] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_lower
	 && trans == blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_upper
					&& trans != blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_lower
	 && trans != blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_upper
					&& trans == blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if ((order == blas_rowmajor && uplo == blas_upper
	 && trans != blas_no_trans) || (order == blas_colmajor
					&& uplo == blas_lower
					&& trans == blas_no_trans)) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  }
}
