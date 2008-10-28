#include "blas_extended.h"
#include "blas_extended_test.h"




void sge_sum_mv_commit(enum blas_order_type order,
		       int m, int n, float *A, int lda, float *y, int row)
{
  int i;
  float *A_i = A;
  float *y_i = y;
  float ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  m_i = m;
  n_i = n;
  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp = y_i[i * inc];
      A_i[(i + lda * row) * inc] = ytemp;
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp = y_i[i * inc];
      A_i[(row + lda * i) * inc] = ytemp;
    }
  }
}
void sge_sum_mv_copy(enum blas_order_type order,
		     int m, int n, float *A, int lda, float *y, int row)
{
  int i;
  float *A_i = A;
  float *y_i = y;
  float ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  m_i = m;
  n_i = n;

  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp = A_i[(i + lda * row) * inc];
      y_i[i * inc] = ytemp;
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp = A_i[(row + lda * i) * inc];
      y_i[i * inc] = ytemp;
    }
  }
}
void dge_sum_mv_commit(enum blas_order_type order,
		       int m, int n, double *A, int lda, double *y, int row)
{
  int i;
  double *A_i = A;
  double *y_i = y;
  double ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  m_i = m;
  n_i = n;
  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp = y_i[i * inc];
      A_i[(i + lda * row) * inc] = ytemp;
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp = y_i[i * inc];
      A_i[(row + lda * i) * inc] = ytemp;
    }
  }
}
void dge_sum_mv_copy(enum blas_order_type order,
		     int m, int n, double *A, int lda, double *y, int row)
{
  int i;
  double *A_i = A;
  double *y_i = y;
  double ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  m_i = m;
  n_i = n;

  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp = A_i[(i + lda * row) * inc];
      y_i[i * inc] = ytemp;
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp = A_i[(row + lda * i) * inc];
      y_i[i * inc] = ytemp;
    }
  }
}
void cge_sum_mv_commit(enum blas_order_type order,
		       int m, int n, void *A, int lda, void *y, int row)
{
  int i;
  float *A_i = (float *) A;
  float *y_i = (float *) y;
  float ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  m_i = m;
  n_i = n;
  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp[0] = y_i[i * inc];
      ytemp[1] = y_i[i * inc + 1];
      A_i[(i + lda * row) * inc] = ytemp[0];
      A_i[(i + lda * row) * inc + 1] = ytemp[1];
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp[0] = y_i[i * inc];
      ytemp[1] = y_i[i * inc + 1];
      A_i[(row + lda * i) * inc] = ytemp[0];
      A_i[(row + lda * i) * inc + 1] = ytemp[1];
    }
  }
}
void cge_sum_mv_copy(enum blas_order_type order,
		     int m, int n, void *A, int lda, void *y, int row)
{
  int i;
  float *A_i = (float *) A;
  float *y_i = (float *) y;
  float ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  m_i = m;
  n_i = n;

  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp[0] = A_i[(i + lda * row) * inc];
      ytemp[1] = A_i[(i + lda * row) * inc + 1];
      y_i[i * inc] = ytemp[0];
      y_i[i * inc + 1] = ytemp[1];
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp[0] = A_i[(row + lda * i) * inc];
      ytemp[1] = A_i[(row + lda * i) * inc + 1];
      y_i[i * inc] = ytemp[0];
      y_i[i * inc + 1] = ytemp[1];
    }
  }
}
void zge_sum_mv_commit(enum blas_order_type order,
		       int m, int n, void *A, int lda, void *y, int row)
{
  int i;
  double *A_i = (double *) A;
  double *y_i = (double *) y;
  double ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  m_i = m;
  n_i = n;
  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp[0] = y_i[i * inc];
      ytemp[1] = y_i[i * inc + 1];
      A_i[(i + lda * row) * inc] = ytemp[0];
      A_i[(i + lda * row) * inc + 1] = ytemp[1];
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp[0] = y_i[i * inc];
      ytemp[1] = y_i[i * inc + 1];
      A_i[(row + lda * i) * inc] = ytemp[0];
      A_i[(row + lda * i) * inc + 1] = ytemp[1];
    }
  }
}
void zge_sum_mv_copy(enum blas_order_type order,
		     int m, int n, void *A, int lda, void *y, int row)
{
  int i;
  double *A_i = (double *) A;
  double *y_i = (double *) y;
  double ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  m_i = m;
  n_i = n;

  if (order == blas_rowmajor) {
    for (i = 0; i < n_i; i++) {
      ytemp[0] = A_i[(i + lda * row) * inc];
      ytemp[1] = A_i[(i + lda * row) * inc + 1];
      y_i[i * inc] = ytemp[0];
      y_i[i * inc + 1] = ytemp[1];
    }
  } else {

    for (i = 0; i < n_i; i++) {
      ytemp[0] = A_i[(row + lda * i) * inc];
      ytemp[1] = A_i[(row + lda * i) * inc + 1];
      y_i[i * inc] = ytemp[0];
      y_i[i * inc + 1] = ytemp[1];
    }
  }
}
