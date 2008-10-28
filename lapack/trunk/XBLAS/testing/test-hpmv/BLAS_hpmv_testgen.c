
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void BLAS_chpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_chpmv{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_chemv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    chpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_chpmv_testgen */
void BLAS_zhpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhpmv{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zhemv_testgen(norm, order, uplo,
		       n, randomize, alpha, alpha_flag, beta, beta_flag,
		       a_full, n /* lda */ , x, incx, y, incy,
		       seed, r_true_l, r_true_t);
    zhpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_zhpmv_testgen */
void BLAS_zhpmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhpmv_c_z{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zhemv_c_z_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    chpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_zhpmv_c_z_testgen */
void BLAS_zhpmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhpmv_z_c{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zhemv_z_c_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    zhpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_zhpmv_z_c_testgen */
void BLAS_zhpmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhpmv_c_c{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) void*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zhemv_c_c_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    chpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_zhpmv_c_c_testgen */
void BLAS_zhpmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_zhpmv_z_d{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) double*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    double *a_full;
    a_full = (double *) blas_malloc(n * n * sizeof(double) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_zhemv_z_d_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    zhpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_zhpmv_z_d_testgen */
void BLAS_chpmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to BLAS_chpmv_c_s{_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) void* The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) float*
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) void*
 *         generated vector y that will be used as an input to HPMV.
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
       R1 = alpha * A1 * x + beta * y1
       R2 = alpha * A2 * x + beta * y2
       where all the matrices and vectors are real.  Then let R = R1 + i R2, 
       A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
       symmetric, and A2 is a skew matrix (trans(A2) = -A2).

       hemv_testgen_body does all that, so just call it, then
       finally, at the end, pack it all properly.
     */

    float *a_full;
    a_full = (float *) blas_malloc(n * n * sizeof(float) * 2);
    if (n * n > 0 && a_full == NULL) {
      BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
    }

    BLAS_chemv_c_s_testgen(norm, order, uplo,
			   n, randomize, alpha, alpha_flag, beta, beta_flag,
			   a_full, n /* lda */ , x, incx, y, incy,
			   seed, r_true_l, r_true_t);
    chpmv_pack_matrix(order, uplo, n, a, a_full, n /*lda */ );

    blas_free(a_full);
  }

}				/* end BLAS_chpmv_c_s_testgen */
