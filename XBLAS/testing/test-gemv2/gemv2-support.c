#include "blas_extended.h"
#include "blas_extended_test.h"






void sgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   float *A, int lda, float *y, int row)
{
  int i;
  float *A_i = A;
  float *y_i = y;
  float ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }
  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void sgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 float *A, int lda, float *y, int row)
{
  int i;
  float *A_i = A;
  float *y_i = y;
  float ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void dgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   double *A, int lda, double *y, int row)
{
  int i;
  double *A_i = A;
  double *y_i = y;
  double ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }
  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void dgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 double *A, int lda, double *y, int row)
{
  int i;
  double *A_i = A;
  double *y_i = y;
  double ytemp;
  int m_i, n_i;
  int inc;

  inc = 1;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void cgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   void *A, int lda, void *y, int row)
{
  int i;
  float *A_i = (float *) A;
  float *y_i = (float *) y;
  float ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }
  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void cgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 void *A, int lda, void *y, int row)
{
  int i;
  float *A_i = (float *) A;
  float *y_i = (float *) y;
  float ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void zgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   void *A, int lda, void *y, int row)
{
  int i;
  double *A_i = (double *) A;
  double *y_i = (double *) y;
  double ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }
  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
void zgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 void *A, int lda, void *y, int row)
{
  int i;
  double *A_i = (double *) A;
  double *y_i = (double *) y;
  double ytemp[2];
  int m_i, n_i;
  int inc;

  inc = 1;
  inc *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  if ((order == blas_rowmajor && trans == blas_no_trans) ||
      (order == blas_colmajor && trans != blas_no_trans)) {
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
