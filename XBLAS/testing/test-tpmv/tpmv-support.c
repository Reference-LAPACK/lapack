#include "blas_extended.h"
#include "blas_extended_test.h"

void stpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const float *a,
		    float *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from a to a_vec
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * a           (input) float*
 *
 * a_vec            (input) float*
 *
 * row          (input) int
 *
 */
{
  int i, j, ind, stride, inc = 1;
  const float *a_i = a;
  float *a_vec_i = a_vec;
  float tmp;
  float tmp2;

  tmp2 = 0.0;
  for (j = 0; j < n; j++) {
    a_vec_i[j * inc] = tmp2;
  }

  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    ind = row * inc;
    stride = (n - 1) * inc;
    for (i = 0; i <= row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += stride;
      stride -= inc;

    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += stride;
      stride += inc;

    }
  }


}

void dtpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const double *a,
		    double *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from a to a_vec
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * a           (input) double*
 *
 * a_vec            (input) double*
 *
 * row          (input) int
 *
 */
{
  int i, j, ind, stride, inc = 1;
  const double *a_i = a;
  double *a_vec_i = a_vec;
  double tmp;
  double tmp2;

  tmp2 = 0.0;
  for (j = 0; j < n; j++) {
    a_vec_i[j * inc] = tmp2;
  }

  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    ind = row * inc;
    stride = (n - 1) * inc;
    for (i = 0; i <= row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += stride;
      stride -= inc;

    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += stride;
      stride += inc;

    }
  }


}

void ctpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const void *a,
		    void *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from a to a_vec
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * a           (input) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, j, ind, stride, inc = 1;
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;
  float tmp[2];
  float tmp2[2];
  inc *= 2;
  tmp2[0] = tmp2[1] = 0.0;
  for (j = 0; j < n; j++) {
    a_vec_i[j * inc] = tmp2[0];
    a_vec_i[j * inc + 1] = tmp2[1];
  }

  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    ind = row * inc;
    stride = (n - 1) * inc;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += stride;
      stride -= inc;

    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += stride;
      stride += inc;

    }
  }


}

void ztpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const void *a,
		    void *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from a to a_vec
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * a           (input) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, j, ind, stride, inc = 1;
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;
  double tmp[2];
  double tmp2[2];
  inc *= 2;
  tmp2[0] = tmp2[1] = 0.0;
  for (j = 0; j < n; j++) {
    a_vec_i[j * inc] = tmp2[0];
    a_vec_i[j * inc + 1] = tmp2[1];
  }

  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    ind = row * inc;
    stride = (n - 1) * inc;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += stride;
      stride -= inc;

    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += stride;
      stride += inc;

    }
  }


}


void stpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, float *a,
		      const float *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a_vec to a
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of a and the length of vector a_vec
 *
 * a            (output) float*
 *
 * a_vec            (input) float*
 *
 * row          (input) int
 *
 */
{
  int i, stride, ind, inc = 1;
  float *a_i = a;
  const float *a_vec_i = a_vec;
  float tmp;





  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    stride = (n - 1) * inc;
    ind = row * inc;
    for (i = 0; i <= row; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += stride;
      stride -= inc;


    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += stride;
      stride += inc;

    }
  }
}

void dtpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, double *a,
		      const double *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a_vec to a
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of a and the length of vector a_vec
 *
 * a            (output) double*
 *
 * a_vec            (input) double*
 *
 * row          (input) int
 *
 */
{
  int i, stride, ind, inc = 1;
  double *a_i = a;
  const double *a_vec_i = a_vec;
  double tmp;





  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    stride = (n - 1) * inc;
    ind = row * inc;
    for (i = 0; i <= row; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += stride;
      stride -= inc;


    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += stride;
      stride += inc;

    }
  }
}

void ctpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, void *a,
		      const void *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a_vec to a
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of a and the length of vector a_vec
 *
 * a            (output) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, stride, ind, inc = 1;
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;
  float tmp[2];

  inc *= 2;



  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    stride = (n - 1) * inc;
    ind = row * inc;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += stride;
      stride -= inc;


    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += stride;
      stride += inc;

    }
  }
}

void ztpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, void *a,
		      const void *a_vec, int row)
/*
 * Purpose
 * =======
 *
 * Copy a_vec to a
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of a and the length of vector a_vec
 *
 * a            (output) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, stride, ind, inc = 1;
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;
  double tmp[2];

  inc *= 2;



  if (((order == blas_rowmajor) && (uplo == blas_upper)
       && (trans != blas_no_trans)) || ((order == blas_colmajor)
					&& (uplo == blas_lower)
					&& (trans == blas_no_trans))) {
    /* colmajor/lower. */
    stride = (n - 1) * inc;
    ind = row * inc;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += stride;
      stride -= inc;


    }

  } else
    if (((order == blas_rowmajor) && (uplo == blas_lower)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_upper)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/lower. */

    ind = row * (row + 1) * inc / 2;
    for (i = 0; i <= row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;

    }


  } else
    if (((order == blas_rowmajor) && (uplo == blas_upper)
	 && (trans == blas_no_trans)) || ((order == blas_colmajor)
					  && (uplo == blas_lower)
					  && (trans != blas_no_trans))) {
    /* Pretend it is rowmajor/upper. */


    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;

    }
  } else {
    /* Pretend it is colmajor/upper. */

    ind = (row + (row * (row + 1)) / 2) * inc;
    stride = (row + 1) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += stride;
      stride += inc;

    }
  }
}
