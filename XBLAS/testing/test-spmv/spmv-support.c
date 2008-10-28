#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"


void sspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const float *a, float *a_vec, int row)
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
 *              Dimension of  and the length of vector x
 *
 * a           (input) float*
 *
 * a_vec            (input) float*
 *
 * row          (input) int
 *
 */
{
  int i, ind, inc = 1;
  const float *a_i = a;
  float *a_vec_i = a_vec;
  float tmp;



  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += (i + 1) * inc;
    }
  }
}


void dspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const double *a, double *a_vec, int row)
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
 *              Dimension of  and the length of vector x
 *
 * a           (input) double*
 *
 * a_vec            (input) double*
 *
 * row          (input) int
 *
 */
{
  int i, ind, inc = 1;
  const double *a_i = a;
  double *a_vec_i = a_vec;
  double tmp;



  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_i[ind];
      a_vec_i[i * inc] = tmp;
      ind += (i + 1) * inc;
    }
  }
}


void cspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *a, void *a_vec, int row)
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
 *              Dimension of  and the length of vector x
 *
 * a           (input) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, ind, inc = 1;
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;
  float tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}


void zspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *a, void *a_vec, int row)
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
 *              Dimension of  and the length of vector x
 *
 * a           (input) void*
 *
 * a_vec            (input) void*
 *
 * row          (input) int
 *
 */
{
  int i, ind, inc = 1;
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;
  double tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}



void sspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, float *a, const float *a_vec, int row)
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
  int i, ind, inc = 1;
  float *a_i = a;
  const float *a_vec_i = a_vec;
  float tmp;



  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is rowmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;
    }
  } else {
    /* Pretend it is colmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += (i + 1) * inc;
    }
  }
}


void dspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double *a, const double *a_vec, int row)
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
  int i, ind, inc = 1;
  double *a_i = a;
  const double *a_vec_i = a_vec;
  double tmp;



  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is rowmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += inc;
    }
  } else {
    /* Pretend it is colmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp = a_vec_i[i * inc];
      a_i[ind] = tmp;
      ind += (i + 1) * inc;
    }
  }
}


void cspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, const void *a_vec, int row)
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
  int i, ind, inc = 1;
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;
  float tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is rowmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is colmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}


void zspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, const void *a_vec, int row)
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
  int i, ind, inc = 1;
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;
  double tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is rowmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is colmajor/upper.  We can do this in the
       symmetric case. */
    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}


void sspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, float *a_packed, float *a_full, int lda)

/*
 *  Packs the he matrix a_full into packed form a.
 */
{
  int row;
  float *a_row;;

  a_row = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && a_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (row = 0; row < n; row++) {
    ssy_copy_row(order, uplo, n, a_full, lda, a_row, row);
    sspmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}
void dspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double *a_packed, double *a_full, int lda)

/*
 *  Packs the he matrix a_full into packed form a.
 */
{
  int row;
  double *a_row;;

  a_row = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && a_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (row = 0; row < n; row++) {
    dsy_copy_row(order, uplo, n, a_full, lda, a_row, row);
    dspmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}
void cspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda)

/*
 *  Packs the he matrix a_full into packed form a.
 */
{
  int row;
  float *a_row;;

  a_row = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (row = 0; row < n; row++) {
    csy_copy_row(order, uplo, n, a_full, lda, a_row, row);
    cspmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}
void zspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda)

/*
 *  Packs the he matrix a_full into packed form a.
 */
{
  int row;
  double *a_row;;

  a_row = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  for (row = 0; row < n; row++) {
    zsy_copy_row(order, uplo, n, a_full, lda, a_row, row);
    zspmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}

void sprint_spmv_matrix(float *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  {
    int row;
    float *x;
    x = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && x == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    for (row = 0; row < n; row++) {
      sspmv_copy_row(order, uplo, n, a, x, row);
      sprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
void dprint_spmv_matrix(double *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  {
    int row;
    double *x;
    x = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && x == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    for (row = 0; row < n; row++) {
      dspmv_copy_row(order, uplo, n, a, x, row);
      dprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
void cprint_spmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  {
    int row;
    float *x;
    x = (float *) blas_malloc(n * sizeof(float) * 2);
    if (n > 0 && x == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    for (row = 0; row < n; row++) {
      cspmv_copy_row(order, uplo, n, a, x, row);
      cprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
void zprint_spmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  {
    int row;
    double *x;
    x = (double *) blas_malloc(n * sizeof(double) * 2);
    if (n > 0 && x == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    for (row = 0; row < n; row++) {
      zspmv_copy_row(order, uplo, n, a, x, row);
      zprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
