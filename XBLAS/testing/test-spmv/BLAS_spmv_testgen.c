#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void BLAS_sspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			float *alpha, int alpha_flag, float *beta,
			int beta_flag, float *a, float *x, int incx, float *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_sspmv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) float*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) float*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) float*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_ssymv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    sspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_sspmv_testgen */
void BLAS_dspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			double *alpha, int alpha_flag, double *beta,
			int beta_flag, double *a, double *x, int incx,
			double *y, int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dspmv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_dsymv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    dspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_dspmv_testgen */
void BLAS_cspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cspmv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_csymv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    cspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_cspmv_testgen */
void BLAS_zspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    zspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_testgen */
void BLAS_cspmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cspmv_s_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_csymv_s_s_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    sspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_cspmv_s_s_testgen */
void BLAS_cspmv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cspmv_s_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_csymv_s_c_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    sspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_cspmv_s_c_testgen */
void BLAS_cspmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_cspmv_c_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_csymv_c_s_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    cspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_cspmv_c_s_testgen */
void BLAS_zspmv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_d_d{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_d_d_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    dspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_d_d_testgen */
void BLAS_zspmv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_d_z{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_d_z_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    dspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_d_z_testgen */
void BLAS_zspmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_z_d{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_z_d_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    zspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_z_d_testgen */
void BLAS_dspmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, float *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dspmv_s_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_dsymv_s_s_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    sspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_dspmv_s_s_testgen */
void BLAS_dspmv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, double *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dspmv_s_d{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) float* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_dsymv_s_d_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    sspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_dspmv_s_d_testgen */
void BLAS_dspmv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, float *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_dspmv_d_s{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) double*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) double*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) double* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) double*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double));
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_dsymv_d_s_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    dspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_dspmv_d_s_testgen */
void BLAS_zspmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_c_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_c_c_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    cspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_c_c_testgen */
void BLAS_zspmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_c_z{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_c_z_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    cspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_c_z_testgen */
void BLAS_zspmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zspmv_z_c{_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) void*
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) void*
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) void* The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */
{
  {

    /* Strategy:  

       Use the SYMV generator, then simply pack the 
       output.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zsymv_z_c_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    zspmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }
}				/* end BLAS_zspmv_z_c_testgen */
