
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void sskew_commit_row_hbmv(enum blas_order_type order,
			   enum blas_uplo_type uplo, int n, float *a, int k,
			   int lda, float *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  float a_elem;
  float *a_i = a;
  const float *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }






  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (uplo == blas_upper) {
      a_elem = -a_elem;;
    }
    a_i[ai] = a_elem;
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (uplo == blas_lower) {
      a_elem = -a_elem;;
    }
    a_i[ai] = a_elem;
  }
}
void dskew_commit_row_hbmv(enum blas_order_type order,
			   enum blas_uplo_type uplo, int n, double *a, int k,
			   int lda, double *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  double a_elem;
  double *a_i = a;
  const double *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }






  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (uplo == blas_upper) {
      a_elem = -a_elem;;
    }
    a_i[ai] = a_elem;
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (uplo == blas_lower) {
      a_elem = -a_elem;;
    }
    a_i[ai] = a_elem;
  }
}

void chbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  float a_elem[2];
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;

  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (uplo == blas_upper) {
      a_elem[1] = -a_elem[1];;
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (uplo == blas_lower) {
      a_elem[1] = -a_elem[1];;
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }
}
void zhbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  double a_elem[2];
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;

  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (uplo == blas_upper) {
      a_elem[1] = -a_elem[1];;
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (uplo == blas_lower) {
      a_elem[1] = -a_elem[1];;
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }
}

void sskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo,
			 int n, float *a, int k, int lda,
			 float *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  float a_elem;
  const float *a_i = a;
  float *a_vec_i = a_vec;


  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }






  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    if (uplo == blas_upper) {
      a_elem = -a_elem;;
    }
    a_vec_i[vi] = a_elem;
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    if (i == row) {

    }
    if (uplo == blas_lower) {
      a_elem = -a_elem;;
    }
    a_vec_i[vi] = a_elem;
  }
  for (; i < n; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
  }


}
void dskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo,
			 int n, double *a, int k, int lda,
			 double *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  double a_elem;
  const double *a_i = a;
  double *a_vec_i = a_vec;


  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }






  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    if (uplo == blas_upper) {
      a_elem = -a_elem;;
    }
    a_vec_i[vi] = a_elem;
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    if (i == row) {

    }
    if (uplo == blas_lower) {
      a_elem = -a_elem;;
    }
    a_vec_i[vi] = a_elem;
  }
  for (; i < n; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
  }


}

void chbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  float a_elem[2];
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;


  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
    a_vec_i[vi + 1] = 0.0;
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (uplo == blas_upper) {
      a_elem[1] = -a_elem[1];;
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (i == row) {
      a_elem[1] = 0.0;
    }
    if (uplo == blas_lower) {
      a_elem[1] = -a_elem[1];;
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
  for (; i < n; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
    a_vec_i[vi + 1] = 0.0;
  }
}
void zhbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  double a_elem[2];
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;


  if (order == blas_colmajor) {
    if (uplo == blas_upper) {
      incai1 = 1;
      incai2 = lda - 1;
    } else {
      incai1 = lda - 1;
      incai2 = 1;
    }
  } else {
    if (uplo == blas_upper) {
      incai1 = lda - 1;
      incai2 = 1;
    } else {
      incai1 = 1;
      incai2 = lda - 1;
    }
  }

  ai = 0;
  if ((uplo == blas_upper && order == blas_colmajor) ||
      (uplo == blas_lower && order == blas_rowmajor)) {
    /* starting place */
    ai = lda * row + ((row < k) ? (k - row) : 0);
  } else {
    /* starting place */
    ai = (row > k) ? (k + lda * (row - k)) : row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row - k; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
    a_vec_i[vi + 1] = 0.0;
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (uplo == blas_upper) {
      a_elem[1] = -a_elem[1];;
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }

  loopmax = MIN(row + k + 1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (i == row) {
      a_elem[1] = 0.0;
    }
    if (uplo == blas_lower) {
      a_elem[1] = -a_elem[1];;
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
  for (; i < n; i++, vi += incvi) {
    a_vec_i[vi] = 0.0;
    a_vec_i[vi + 1] = 0.0;
  }
}

void cprint_hbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  int row;
  float *x;
  x = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (row = 0; row < n; row++) {
    chbmv_copy_row(order, uplo, n, a, k, lda, x, row);
    cprint_vector(x, n, 1, NULL);
  }
  printf("\n");
  blas_free(x);

}
void zprint_hbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{

  int row;
  double *x;
  x = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (row = 0; row < n; row++) {
    zhbmv_copy_row(order, uplo, n, a, k, lda, x, row);
    zprint_vector(x, n, 1, NULL);
  }
  printf("\n");
  blas_free(x);

}
