
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


void sgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, float *AB, int lda, float *y,
		   int row, int *nfix2, int *nmix, int *ysize)
{
  int ra, la, i, iy;
  int n_i;
  int inc = 1;
  float *y_i = y;
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else if ((order == blas_colmajor) && (trans == blas_trans)) {
    ra = kl;
    la = ku;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else {			/* rowmajor and blas_trans */
    ra = kl;
    la = ku;
  }

  if (trans == blas_no_trans)
    n_i = n;
  else
    n_i = m;

  /* DETERMINE SIZE OF VECTOR TO BE SUBMITTED TO TESTGEN */

  if (row + ra + 1 < n_i) {
    *ysize = ra + row + 1;
  } else {
    *ysize = n_i;
  }

  /* SET NMIX AND NFIX */
  if (row == 0) {
    *nmix = 0;
    *nfix2 = 0;
  } else {
    if (row > la) {
      *nfix2 = row - la;
    } else {
      *nfix2 = 0;
    }
    if (*nfix2 > n_i)
      *nfix2 = n_i;
    if (row <= n_i - ra - 1) {
      *nmix = (*ysize) - (*nfix2) - 1;
    } else {
      *nmix = (*ysize) - (*nfix2);
    }
  }

  /* SET ALL VALUES OF Y TO  = 0.0; */

  for (i = 0, iy = 0; i < *ysize; i++, iy += inc) {
    y_i[iy] = 0.0;
  }
}
void sgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, float *AB, int lda, float *y,
		  int row)
{
  int i, j;
  float *AB_i = AB;
  float *y_i = y;
  float ytemp;
  int inc = 1;

  row++;


  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = y_i[(j - 1) * inc];
	AB_i[(ku + row - j + lda * (j - 1)) * inc] = ytemp;
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = y_i[(j - 1) * inc];
	AB_i[((row - 1) * lda + kl + j - row) * inc] = ytemp;
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = y_i[(i - 1) * inc];
	AB_i[(ku + i - row + lda * (row - 1)) * inc] = ytemp;
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = y_i[(i - 1) * inc];
	AB_i[((i - 1) * lda + kl - i + row) * inc] = ytemp;
      }
    }
  }
}
void sgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const float *AB, int lda, float *y,
		int row)
{
  int i, j, iy;
  int max_mn;
  const float *AB_i = AB;
  float *y_i = y;
  float ytemp;
  int inc = 1;

  row++;

  max_mn = MAX(m, n);

  for (i = 0, iy = 0; i < max_mn; i++, iy += inc) {
    y_i[iy] = 0.0;
  }

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = AB_i[(ku + row - j + lda * (j - 1)) * inc];
	y_i[(j - 1) * inc] = ytemp;
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = AB_i[((row - 1) * lda + kl + j - row) * inc];
	y_i[(j - 1) * inc] = ytemp;
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = AB_i[((ku + i - row) + lda * (row - 1)) * inc];
	y_i[(i - 1) * inc] = ytemp;
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = AB_i[((i - 1) * lda + kl - i + row) * inc];
	y_i[(i - 1) * inc] = ytemp;
      }
    }
  }
}

void cgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, void *AB, int lda, void *y,
		   int row, int *nfix2, int *nmix, int *ysize)
{
  int ra, la, i, iy;
  int n_i;
  int inc = 1;
  float *y_i = (float *) y;
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else if ((order == blas_colmajor) && (trans == blas_trans)) {
    ra = kl;
    la = ku;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else {			/* rowmajor and blas_trans */
    ra = kl;
    la = ku;
  }

  if (trans == blas_no_trans)
    n_i = n;
  else
    n_i = m;

  /* DETERMINE SIZE OF VECTOR TO BE SUBMITTED TO TESTGEN */

  if (row + ra + 1 < n_i) {
    *ysize = ra + row + 1;
  } else {
    *ysize = n_i;
  }

  /* SET NMIX AND NFIX */
  if (row == 0) {
    *nmix = 0;
    *nfix2 = 0;
  } else {
    if (row > la) {
      *nfix2 = row - la;
    } else {
      *nfix2 = 0;
    }
    if (*nfix2 > n_i)
      *nfix2 = n_i;
    if (row <= n_i - ra - 1) {
      *nmix = (*ysize) - (*nfix2) - 1;
    } else {
      *nmix = (*ysize) - (*nfix2);
    }
  }

  /* SET ALL VALUES OF Y TO  = 0.0; */
  inc *= 2;
  for (i = 0, iy = 0; i < *ysize; i++, iy += inc) {
    y_i[iy] = 0.0;
    y_i[iy + 1] = 0.0;
  }
}
void cgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, void *AB, int lda, void *y,
		  int row)
{
  int i, j;
  float *AB_i = (float *) AB;
  float *y_i = (float *) y;
  float ytemp[2];
  int inc = 1;

  row++;
  inc *= 2;

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = y_i[(j - 1) * inc];
	ytemp[1] = y_i[(j - 1) * inc + 1];
	AB_i[(ku + row - j + lda * (j - 1)) * inc] = ytemp[0];
	AB_i[(ku + row - j + lda * (j - 1)) * inc + 1] = ytemp[1];
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = y_i[(j - 1) * inc];
	ytemp[1] = y_i[(j - 1) * inc + 1];
	AB_i[((row - 1) * lda + kl + j - row) * inc] = ytemp[0];
	AB_i[((row - 1) * lda + kl + j - row) * inc + 1] = ytemp[1];
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = y_i[(i - 1) * inc];
	ytemp[1] = y_i[(i - 1) * inc + 1];
	AB_i[(ku + i - row + lda * (row - 1)) * inc] = ytemp[0];
	AB_i[(ku + i - row + lda * (row - 1)) * inc + 1] = ytemp[1];
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = y_i[(i - 1) * inc];
	ytemp[1] = y_i[(i - 1) * inc + 1];
	AB_i[((i - 1) * lda + kl - i + row) * inc] = ytemp[0];
	AB_i[((i - 1) * lda + kl - i + row) * inc + 1] = ytemp[1];
      }
    }
  }
}
void cgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const void *AB, int lda, void *y,
		int row)
{
  int i, j, iy;
  int max_mn;
  const float *AB_i = (float *) AB;
  float *y_i = (float *) y;
  float ytemp[2];
  int inc = 1;

  row++;
  inc *= 2;
  max_mn = MAX(m, n);

  for (i = 0, iy = 0; i < max_mn; i++, iy += inc) {
    y_i[iy] = 0.0;
    y_i[iy + 1] = 0.0;
  }

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = AB_i[(ku + row - j + lda * (j - 1)) * inc];
	ytemp[1] = AB_i[(ku + row - j + lda * (j - 1)) * inc + 1];
	y_i[(j - 1) * inc] = ytemp[0];
	y_i[(j - 1) * inc + 1] = ytemp[1];
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = AB_i[((row - 1) * lda + kl + j - row) * inc];
	ytemp[1] = AB_i[((row - 1) * lda + kl + j - row) * inc + 1];
	y_i[(j - 1) * inc] = ytemp[0];
	y_i[(j - 1) * inc + 1] = ytemp[1];
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = AB_i[((ku + i - row) + lda * (row - 1)) * inc];
	ytemp[1] = AB_i[((ku + i - row) + lda * (row - 1)) * inc + 1];
	y_i[(i - 1) * inc] = ytemp[0];
	y_i[(i - 1) * inc + 1] = ytemp[1];
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = AB_i[((i - 1) * lda + kl - i + row) * inc];
	ytemp[1] = AB_i[((i - 1) * lda + kl - i + row) * inc + 1];
	y_i[(i - 1) * inc] = ytemp[0];
	y_i[(i - 1) * inc + 1] = ytemp[1];
      }
    }
  }
}

void dgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, double *AB, int lda,
		   double *y, int row, int *nfix2, int *nmix, int *ysize)
{
  int ra, la, i, iy;
  int n_i;
  int inc = 1;
  double *y_i = y;
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else if ((order == blas_colmajor) && (trans == blas_trans)) {
    ra = kl;
    la = ku;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else {			/* rowmajor and blas_trans */
    ra = kl;
    la = ku;
  }

  if (trans == blas_no_trans)
    n_i = n;
  else
    n_i = m;

  /* DETERMINE SIZE OF VECTOR TO BE SUBMITTED TO TESTGEN */

  if (row + ra + 1 < n_i) {
    *ysize = ra + row + 1;
  } else {
    *ysize = n_i;
  }

  /* SET NMIX AND NFIX */
  if (row == 0) {
    *nmix = 0;
    *nfix2 = 0;
  } else {
    if (row > la) {
      *nfix2 = row - la;
    } else {
      *nfix2 = 0;
    }
    if (*nfix2 > n_i)
      *nfix2 = n_i;
    if (row <= n_i - ra - 1) {
      *nmix = (*ysize) - (*nfix2) - 1;
    } else {
      *nmix = (*ysize) - (*nfix2);
    }
  }

  /* SET ALL VALUES OF Y TO  = 0.0; */

  for (i = 0, iy = 0; i < *ysize; i++, iy += inc) {
    y_i[iy] = 0.0;
  }
}
void dgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, double *AB, int lda,
		  double *y, int row)
{
  int i, j;
  double *AB_i = AB;
  double *y_i = y;
  double ytemp;
  int inc = 1;

  row++;


  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = y_i[(j - 1) * inc];
	AB_i[(ku + row - j + lda * (j - 1)) * inc] = ytemp;
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = y_i[(j - 1) * inc];
	AB_i[((row - 1) * lda + kl + j - row) * inc] = ytemp;
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = y_i[(i - 1) * inc];
	AB_i[(ku + i - row + lda * (row - 1)) * inc] = ytemp;
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = y_i[(i - 1) * inc];
	AB_i[((i - 1) * lda + kl - i + row) * inc] = ytemp;
      }
    }
  }
}
void dgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const double *AB, int lda, double *y,
		int row)
{
  int i, j, iy;
  int max_mn;
  const double *AB_i = AB;
  double *y_i = y;
  double ytemp;
  int inc = 1;

  row++;

  max_mn = MAX(m, n);

  for (i = 0, iy = 0; i < max_mn; i++, iy += inc) {
    y_i[iy] = 0.0;
  }

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = AB_i[(ku + row - j + lda * (j - 1)) * inc];
	y_i[(j - 1) * inc] = ytemp;
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp = AB_i[((row - 1) * lda + kl + j - row) * inc];
	y_i[(j - 1) * inc] = ytemp;
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = AB_i[((ku + i - row) + lda * (row - 1)) * inc];
	y_i[(i - 1) * inc] = ytemp;
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp = AB_i[((i - 1) * lda + kl - i + row) * inc];
	y_i[(i - 1) * inc] = ytemp;
      }
    }
  }
}

void zgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, void *AB, int lda, void *y,
		   int row, int *nfix2, int *nmix, int *ysize)
{
  int ra, la, i, iy;
  int n_i;
  int inc = 1;
  double *y_i = (double *) y;
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else if ((order == blas_colmajor) && (trans == blas_trans)) {
    ra = kl;
    la = ku;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else {			/* rowmajor and blas_trans */
    ra = kl;
    la = ku;
  }

  if (trans == blas_no_trans)
    n_i = n;
  else
    n_i = m;

  /* DETERMINE SIZE OF VECTOR TO BE SUBMITTED TO TESTGEN */

  if (row + ra + 1 < n_i) {
    *ysize = ra + row + 1;
  } else {
    *ysize = n_i;
  }

  /* SET NMIX AND NFIX */
  if (row == 0) {
    *nmix = 0;
    *nfix2 = 0;
  } else {
    if (row > la) {
      *nfix2 = row - la;
    } else {
      *nfix2 = 0;
    }
    if (*nfix2 > n_i)
      *nfix2 = n_i;
    if (row <= n_i - ra - 1) {
      *nmix = (*ysize) - (*nfix2) - 1;
    } else {
      *nmix = (*ysize) - (*nfix2);
    }
  }

  /* SET ALL VALUES OF Y TO  = 0.0; */
  inc *= 2;
  for (i = 0, iy = 0; i < *ysize; i++, iy += inc) {
    y_i[iy] = 0.0;
    y_i[iy + 1] = 0.0;
  }
}
void zgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, void *AB, int lda, void *y,
		  int row)
{
  int i, j;
  double *AB_i = (double *) AB;
  double *y_i = (double *) y;
  double ytemp[2];
  int inc = 1;

  row++;
  inc *= 2;

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = y_i[(j - 1) * inc];
	ytemp[1] = y_i[(j - 1) * inc + 1];
	AB_i[(ku + row - j + lda * (j - 1)) * inc] = ytemp[0];
	AB_i[(ku + row - j + lda * (j - 1)) * inc + 1] = ytemp[1];
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = y_i[(j - 1) * inc];
	ytemp[1] = y_i[(j - 1) * inc + 1];
	AB_i[((row - 1) * lda + kl + j - row) * inc] = ytemp[0];
	AB_i[((row - 1) * lda + kl + j - row) * inc + 1] = ytemp[1];
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = y_i[(i - 1) * inc];
	ytemp[1] = y_i[(i - 1) * inc + 1];
	AB_i[(ku + i - row + lda * (row - 1)) * inc] = ytemp[0];
	AB_i[(ku + i - row + lda * (row - 1)) * inc + 1] = ytemp[1];
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = y_i[(i - 1) * inc];
	ytemp[1] = y_i[(i - 1) * inc + 1];
	AB_i[((i - 1) * lda + kl - i + row) * inc] = ytemp[0];
	AB_i[((i - 1) * lda + kl - i + row) * inc + 1] = ytemp[1];
      }
    }
  }
}
void zgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const void *AB, int lda, void *y,
		int row)
{
  int i, j, iy;
  int max_mn;
  const double *AB_i = (double *) AB;
  double *y_i = (double *) y;
  double ytemp[2];
  int inc = 1;

  row++;
  inc *= 2;
  max_mn = MAX(m, n);

  for (i = 0, iy = 0; i < max_mn; i++, iy += inc) {
    y_i[iy] = 0.0;
    y_i[iy + 1] = 0.0;
  }

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = AB_i[(ku + row - j + lda * (j - 1)) * inc];
	ytemp[1] = AB_i[(ku + row - j + lda * (j - 1)) * inc + 1];
	y_i[(j - 1) * inc] = ytemp[0];
	y_i[(j - 1) * inc + 1] = ytemp[1];
      }
    } else {
      for (j = MAX(1, row - kl); j <= MIN(n, row + ku); j++) {
	ytemp[0] = AB_i[((row - 1) * lda + kl + j - row) * inc];
	ytemp[1] = AB_i[((row - 1) * lda + kl + j - row) * inc + 1];
	y_i[(j - 1) * inc] = ytemp[0];
	y_i[(j - 1) * inc + 1] = ytemp[1];
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = AB_i[((ku + i - row) + lda * (row - 1)) * inc];
	ytemp[1] = AB_i[((ku + i - row) + lda * (row - 1)) * inc + 1];
	y_i[(i - 1) * inc] = ytemp[0];
	y_i[(i - 1) * inc + 1] = ytemp[1];
      }
    } else {
      for (i = MAX(1, row - ku); i <= MIN(m, row + kl); i++) {
	ytemp[0] = AB_i[((i - 1) * lda + kl - i + row) * inc];
	ytemp[1] = AB_i[((i - 1) * lda + kl - i + row) * inc + 1];
	y_i[(i - 1) * inc] = ytemp[0];
	y_i[(i - 1) * inc + 1] = ytemp[1];
      }
    }
  }
}
