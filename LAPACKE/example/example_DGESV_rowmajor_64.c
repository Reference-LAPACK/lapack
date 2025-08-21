/*
   LAPACKE_dgesv Example
   =====================

   The program computes the solution to the system of linear
   equations with a square matrix A and multiple
   right-hand sides B, where A is the coefficient matrix
   and b is the right-hand side matrix:

   Description
   ===========

   The routine solves for X the system of linear equations A*X = B,
   where A is an n-by-n matrix, the columns of matrix B are individual
   right-hand sides, and the columns of X are the corresponding
   solutions.

   The LU decomposition with partial pivoting and row interchanges is
   used to factor A as A = P*L*U, where P is a permutation matrix, L
   is unit lower triangular, and U is upper triangular. The factored
   form of A is then used to solve the system of equations A*X = B.

   LAPACKE Interface
   =================

   LAPACKE_dgesv (row-major, high-level) Example Program Results

  -- LAPACKE Example routine --
  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lapacke_64.h>
#include "lapacke_example_aux.h"

/* Main program */
int main(int argc, char **argv) {

        /* Locals */
        int64_t n, nrhs, lda, ldb, info;
        int i, j;
        /* Local arrays */
        double *A, *b;
        int64_t *ipiv;

        /* Default Value */
        n = 5; nrhs = 1;

        /* Arguments */
        for( i = 1; i < argc; i++ ) {
            if( strcmp( argv[i], "-n" ) == 0 ) {
                n  = atoi(argv[i+1]);
                i++;
            }
            if( strcmp( argv[i], "-nrhs" ) == 0 ) {
                nrhs  = atoi(argv[i+1]);
                i++;
            }
        }

        /* Initialization */
        lda=n, ldb=nrhs;
        A = (double *)malloc(n*n*sizeof(double)) ;
        if (A==NULL){ printf("error of memory allocation\n"); exit(0); }
        b = (double *)malloc(n*nrhs*sizeof(double)) ;
        if (b==NULL){ printf("error of memory allocation\n"); free(A); exit(0); }
        ipiv = (int64_t *)malloc(n*sizeof(int64_t)) ;
        if (ipiv==NULL){ printf("error of memory allocation\n"); free(A); free(b); exit(0); }

        for( i = 0; i < n; i++ ) {
                for( j = 0; j < n; j++ ) A[i*lda+j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
        }

        for(i=0;i<n*nrhs;i++)
            b[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;

        /* Print Entry Matrix */
        print_matrix_rowmajor_64( "Entry Matrix A", n, n, A, lda );
        /* Print Right Rand Side */
        print_matrix_rowmajor_64( "Right Rand Side b", n, nrhs, b, ldb );
        printf( "\n" );
        /* Executable statements */
        printf( "LAPACKE_dgesv_64 (row-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        info = LAPACKE_dgesv_64( LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv,
                        b, ldb );
        /* Check for the exact singularity */
        if( info > 0 ) {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%d,%d) is zero, so that A is singular;\n", (int) info, (int) info );
            printf( "the solution could not be computed.\n" );
            free(A);
            free(b);
            free(ipiv);
            exit( 1 );
        }
        if (info <0) {
            free(A);
            free(b);
            free(ipiv);
            exit( 1 );
        }
        /* Print solution */
        print_matrix_rowmajor_64( "Solution", n, nrhs, b, ldb );
        /* Print details of LU factorization */
        print_matrix_rowmajor_64( "Details of LU factorization", n, n, A, lda );
        /* Print pivot indices */
        print_vector_64( "Pivot indices", n, ipiv );
        free(A);
        free(b);
        free(ipiv);
        exit( 0 );
} /* End of LAPACKE_dgesv Example */

