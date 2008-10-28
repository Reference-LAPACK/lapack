#include "blas_extended.h"
#include "blas_extended_test.h"

void str_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const float *T, int ldt,
		  float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
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
 *              Dimension of AP and the length of vector y
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension of matrix T
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  float tmp;
  const float *T_i = T;
  float *y_i = y;


  for (i = 0; i < n * inc; i += inc) {
    y_i[i] = 0.0;
  }

  if ((order == blas_rowmajor &&
       trans == blas_no_trans && uplo == blas_upper) ||
      (order == blas_colmajor &&
       trans != blas_no_trans && uplo == blas_lower)) {
    for (i = row; i < n; i++) {
      tmp = T_i[(row * ldt + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans == blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans != blas_no_trans && uplo == blas_upper)) {
    for (i = 0; i <= row; i++) {
      tmp = T_i[(row * ldt + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_upper)) {
    for (i = row; i < n; i++) {
      tmp = T_i[(i * ldt + row) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_upper) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_lower)) {
    for (i = 0; i <= row; i++) {
      tmp = T_i[(i * ldt + row) * inc];
      y_i[i * inc] = tmp;
    }
  }
}

void dtr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const double *T, int ldt,
		  double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
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
 *              Dimension of AP and the length of vector y
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension of matrix T
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  double tmp;
  const double *T_i = T;
  double *y_i = y;


  for (i = 0; i < n * inc; i += inc) {
    y_i[i] = 0.0;
  }

  if ((order == blas_rowmajor &&
       trans == blas_no_trans && uplo == blas_upper) ||
      (order == blas_colmajor &&
       trans != blas_no_trans && uplo == blas_lower)) {
    for (i = row; i < n; i++) {
      tmp = T_i[(row * ldt + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans == blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans != blas_no_trans && uplo == blas_upper)) {
    for (i = 0; i <= row; i++) {
      tmp = T_i[(row * ldt + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_upper)) {
    for (i = row; i < n; i++) {
      tmp = T_i[(i * ldt + row) * inc];
      y_i[i * inc] = tmp;
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_upper) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_lower)) {
    for (i = 0; i <= row; i++) {
      tmp = T_i[(i * ldt + row) * inc];
      y_i[i * inc] = tmp;
    }
  }
}

void ctr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const void *T, int ldt,
		  void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
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
 *              Dimension of AP and the length of vector y
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension of matrix T
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  float tmp[2];
  const float *T_i = (float *) T;
  float *y_i = (float *) y;
  inc *= 2;

  for (i = 0; i < n * inc; i += inc) {
    y_i[i] = 0.0;
    y_i[i + 1] = 0.0;
  }

  if ((order == blas_rowmajor &&
       trans == blas_no_trans && uplo == blas_upper) ||
      (order == blas_colmajor &&
       trans != blas_no_trans && uplo == blas_lower)) {
    for (i = row; i < n; i++) {
      tmp[0] = T_i[(row * ldt + i) * inc];
      tmp[1] = T_i[(row * ldt + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans == blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans != blas_no_trans && uplo == blas_upper)) {
    for (i = 0; i <= row; i++) {
      tmp[0] = T_i[(row * ldt + i) * inc];
      tmp[1] = T_i[(row * ldt + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_upper)) {
    for (i = row; i < n; i++) {
      tmp[0] = T_i[(i * ldt + row) * inc];
      tmp[1] = T_i[(i * ldt + row) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_upper) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_lower)) {
    for (i = 0; i <= row; i++) {
      tmp[0] = T_i[(i * ldt + row) * inc];
      tmp[1] = T_i[(i * ldt + row) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  }
}

void ztr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const void *T, int ldt,
		  void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
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
 *              Dimension of AP and the length of vector y
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension of matrix T
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  double tmp[2];
  const double *T_i = (double *) T;
  double *y_i = (double *) y;
  inc *= 2;

  for (i = 0; i < n * inc; i += inc) {
    y_i[i] = 0.0;
    y_i[i + 1] = 0.0;
  }

  if ((order == blas_rowmajor &&
       trans == blas_no_trans && uplo == blas_upper) ||
      (order == blas_colmajor &&
       trans != blas_no_trans && uplo == blas_lower)) {
    for (i = row; i < n; i++) {
      tmp[0] = T_i[(row * ldt + i) * inc];
      tmp[1] = T_i[(row * ldt + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans == blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans != blas_no_trans && uplo == blas_upper)) {
    for (i = 0; i <= row; i++) {
      tmp[0] = T_i[(row * ldt + i) * inc];
      tmp[1] = T_i[(row * ldt + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_lower) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_upper)) {
    for (i = row; i < n; i++) {
      tmp[0] = T_i[(i * ldt + row) * inc];
      tmp[1] = T_i[(i * ldt + row) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else if ((order == blas_rowmajor &&
	      trans != blas_no_trans && uplo == blas_upper) ||
	     (order == blas_colmajor &&
	      trans == blas_no_trans && uplo == blas_lower)) {
    for (i = 0; i <= row; i++) {
      tmp[0] = T_i[(i * ldt + row) * inc];
      tmp[1] = T_i[(i * ldt + row) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  }
}
