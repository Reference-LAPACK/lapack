#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void stbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const float *T,
		int ldt, float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 *      NOTE: this function just encapsulates sgbmv_copy.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T. 
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;


  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }
  sgbmv_copy(order, new_trans, n, n, kl, ku, T, ldt, y, row);



}

void stbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, float *T, int ldt,
		  float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 *      NOTE: this function just encapsulates sgbmv_commit.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;


  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }



  sgbmv_commit(order, new_trans, n, n, kl, ku, T, ldt, y, row);




}

void dtbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const double *T,
		int ldt, double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 *      NOTE: this function just encapsulates dgbmv_copy.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T. 
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;


  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }
  dgbmv_copy(order, new_trans, n, n, kl, ku, T, ldt, y, row);



}

void dtbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, double *T,
		  int ldt, double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 *      NOTE: this function just encapsulates dgbmv_commit.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;


  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }



  dgbmv_commit(order, new_trans, n, n, kl, ku, T, ldt, y, row);




}

void ctbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const void *T,
		int ldt, void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 *      NOTE: this function just encapsulates cgbmv_copy.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T. 
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;
  float *y_i = (float *) y;

  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }
  cgbmv_copy(order, new_trans, n, n, kl, ku, T, ldt, y, row);

  if (conj) {
    float y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

}

void ctbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, void *T, int ldt,
		  void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 *      NOTE: this function just encapsulates cgbmv_commit.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;
  float *y_i = (float *) y;

  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }

  if (conj) {
    float y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

  cgbmv_commit(order, new_trans, n, n, kl, ku, T, ldt, y, row);


  if (conj) {
    /* now conjugate back - leave original */
    float y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

}

void ztbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const void *T,
		int ldt, void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 *      NOTE: this function just encapsulates zgbmv_copy.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T. 
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;
  double *y_i = (double *) y;

  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }
  zgbmv_copy(order, new_trans, n, n, kl, ku, T, ldt, y, row);

  if (conj) {
    double y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

}

void ztbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, void *T, int ldt,
		  void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 *      NOTE: this function just encapsulates zgbmv_commit.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  enum blas_trans_type new_trans;
  int conj = 0;
  int kl, ku;
  double *y_i = (double *) y;

  if (uplo == blas_upper) {
    ku = k;
    kl = 0;
  } else {
    ku = 0;
    kl = k;
  }


  if (trans == blas_no_trans) {
    new_trans = blas_no_trans;
  } else if (trans == blas_conj_trans) {
    new_trans = blas_trans;
    conj = 1;
  } else if (trans == blas_trans) {
    new_trans = blas_trans;
  } else {
    /* conj, no trans */
    new_trans = blas_no_trans;
    conj = 1;
  }

  if (conj) {
    double y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

  zgbmv_commit(order, new_trans, n, n, kl, ku, T, ldt, y, row);


  if (conj) {
    /* now conjugate back - leave original */
    double y_elem[2];
    int i, incyi, ni;

    incyi = 1;
    ni = n;
    ni *= 2;
    incyi *= 2;
    for (i = 0; i < ni; i += incyi) {
      y_elem[0] = y_i[i];
      y_elem[1] = y_i[i + 1];
      y_elem[1] = -y_elem[1];
      y_i[i] = y_elem[0];
      y_i[i + 1] = y_elem[1];
    }
  }

}


void sprint_tbsv_matrix(float *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans)
/*
 * Purpose
 * =======
 *
 * Print T
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 */
{
  int i;
  float *T_row;
  T_row = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && T_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (i = 0; i < n; i++) {
    stbsv_copy(order, uplo, trans, n, k, T, ldt, T_row, i);
    sprint_vector(T_row, n, 1, NULL);
  }
  printf("\n");
  blas_free(T_row);
}

void dprint_tbsv_matrix(double *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans)
/*
 * Purpose
 * =======
 *
 * Print T
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 */
{
  int i;
  double *T_row;
  T_row = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && T_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (i = 0; i < n; i++) {
    dtbsv_copy(order, uplo, trans, n, k, T, ldt, T_row, i);
    dprint_vector(T_row, n, 1, NULL);
  }
  printf("\n");
  blas_free(T_row);
}

void cprint_tbsv_matrix(void *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans)
/*
 * Purpose
 * =======
 *
 * Print T
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 */
{
  int i;
  float *T_row;
  T_row = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && T_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (i = 0; i < n; i++) {
    ctbsv_copy(order, uplo, trans, n, k, T, ldt, T_row, i);
    cprint_vector(T_row, n, 1, NULL);
  }
  printf("\n");
  blas_free(T_row);
}

void zprint_tbsv_matrix(void *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans)
/*
 * Purpose
 * =======
 *
 * Print T
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 */
{
  int i;
  double *T_row;
  T_row = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && T_row == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  for (i = 0; i < n; i++) {
    ztbsv_copy(order, uplo, trans, n, k, T, ldt, T_row, i);
    zprint_vector(T_row, n, 1, NULL);
  }
  printf("\n");
  blas_free(T_row);
}
