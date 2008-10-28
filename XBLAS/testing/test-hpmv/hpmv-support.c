#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void chpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, void *a_vec, int row)

/*
 *  Copies the given row of packed hermitianmatrix a into the supplied vector.
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
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      if (i == row) {
	tmp[1] = 0.0;
      }
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
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      if (i == row) {
	tmp[1] = 0.0;;
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}

void zhpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, void *a_vec, int row)

/*
 *  Copies the given row of packed hermitianmatrix a into the supplied vector.
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
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      if (i == row) {
	tmp[1] = 0.0;
      }
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
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_i[ind];
      tmp[1] = a_i[ind + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      if (i == row) {
	tmp[1] = 0.0;;
      }
      a_vec_i[i * inc] = tmp[0];
      a_vec_i[i * inc + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}


void chpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, void *a_vec, int row)

/*
 *  Commits the supplied vector a_vec to the given row of 
 *      packed hermitian matrix a.
 */
{
  int i, ind, inc = 1;
  float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;
  float tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}

void zhpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, void *a_vec, int row)

/*
 *  Commits the supplied vector a_vec to the given row of 
 *      packed hermitian matrix a.
 */
{
  int i, ind, inc = 1;
  double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;
  double tmp[2];

  inc *= 2;

  if (((order == blas_rowmajor) && (uplo == blas_upper)) ||
      ((order == blas_colmajor) && (uplo == blas_lower))) {
    /* Pretend it is colmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * inc;
    for (i = 0; i < row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (n - i - 1) * inc;
    }

    ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }
  } else {
    /* Pretend it is rowmajor/lower.  We can do this in the
       symmetric case. */
    ind = row * (row + 1) * inc / 2;
    for (i = 0; i < row; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_upper) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += inc;
    }

    ind = (row + (row * (row + 1)) / 2) * inc;
    for (i = row; i < n; i++) {
      tmp[0] = a_vec_i[i * inc];
      tmp[1] = a_vec_i[i * inc + 1];
      if (uplo == blas_lower) {
	tmp[1] = -tmp[1];
      }
      a_i[ind] = tmp[0];
      a_i[ind + 1] = tmp[1];
      ind += (i + 1) * inc;
    }
  }
}


void chpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
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
    che_copy_row(order, uplo, blas_left_side, n, a_full, lda, a_row, row);
    chpmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}
void zhpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
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
    zhe_copy_row(order, uplo, blas_left_side, n, a_full, lda, a_row, row);
    zhpmv_commit_row(order, uplo, n, a_packed, a_row, row);
  }

  blas_free(a_row);
}

void cprint_hpmv_matrix(void *a, int n,
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
      chpmv_copy_row(order, uplo, n, a, x, row);
      cprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
void zprint_hpmv_matrix(void *a, int n,
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
      zhpmv_copy_row(order, uplo, n, a, x, row);
      zprint_vector(x, n, 1, NULL);
    }
    printf("\n");
    blas_free(x);
  }

}
