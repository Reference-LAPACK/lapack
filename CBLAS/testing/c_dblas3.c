/*
 *     Written by D.P. Manley, Digital Equipment Corporation.
 *     Prefixed "C_" to BLAS routines and their declarations.
 *
 *     Modified by T. H. Do, 2/19/98, SGI/CRAY Research.
 */
#include <stdlib.h>
#include "cblas.h"
#include "cblas_test.h"
#define  TEST_COL_MJR	0
#define  TEST_ROW_MJR	1
#define  UNDEFINED     -1

void F77_dgemm(CBLAS_INT *layout, char *transpa, char *transpb, CBLAS_INT *m, CBLAS_INT *n,
              CBLAS_INT *k, double *alpha, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb,
              double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN transpa_len, FORTRAN_STRLEN transpb_len
#endif
) {

  double *A, *B, *C;
  CBLAS_INT i,j,LDA, LDB, LDC;
  CBLAS_TRANSPOSE transa, transb;

  get_transpose_type(transpa, &transa);
  get_transpose_type(transpb, &transb);

  if (*layout == TEST_ROW_MJR) {
     if (transa == CblasNoTrans) {
        LDA = *k+1;
        A = (double *)malloc( (*m)*LDA*sizeof( double ) );
        for( i=0; i<*m; i++ )
           for( j=0; j<*k; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else {
        LDA = *m+1;
        A   = ( double* )malloc( LDA*(*k)*sizeof( double ) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*m; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     if (transb == CblasNoTrans) {
        LDB = *n+1;
        B   = ( double* )malloc( (*k)*LDB*sizeof( double ) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ )
              B[i*LDB+j]=b[j*(*ldb)+i];
     }
     else {
        LDB = *k+1;
        B   = ( double* )malloc( LDB*(*n)*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ )
              B[i*LDB+j]=b[j*(*ldb)+i];
     }
     LDC = *n+1;
     C   = ( double* )malloc( (*m)*LDC*sizeof( double ) );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           C[i*LDC+j]=c[j*(*ldc)+i];

     cblas_dgemm( CblasRowMajor, transa, transb, *m, *n, *k, *alpha, A, LDA,
                  B, LDB, *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dgemm( CblasColMajor, transa, transb, *m, *n, *k, *alpha, a, *lda,
                  b, *ldb, *beta, c, *ldc );
  else
     cblas_dgemm( UNDEFINED, transa, transb, *m, *n, *k, *alpha, a, *lda,
                  b, *ldb, *beta, c, *ldc );
}

void F77_dgemmtr(CBLAS_INT *layout, char *uplop, char *transpa, char *transpb, CBLAS_INT *n,
     CBLAS_INT *k, double *alpha, double *a, CBLAS_INT *lda,
     double *b, CBLAS_INT *ldb, double *beta,
     double *c, CBLAS_INT *ldc ) {

  double *A, *B, *C;
  CBLAS_INT i,j,LDA, LDB, LDC;
  CBLAS_TRANSPOSE transa, transb;
  CBLAS_UPLO uplo;

  get_transpose_type(transpa, &transa);
  get_transpose_type(transpb, &transb);
  get_uplo_type(uplop, &uplo);

  if (*layout == TEST_ROW_MJR) {
     if (transa == CblasNoTrans) {
        LDA = *k+1;
        A=(double*)malloc((*n)*LDA*sizeof(double));
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ ) {
              A[i*LDA+j]=a[j*(*lda)+i];
           }
     }
     else {
        LDA = *n+1;
        A=(double* )malloc(LDA*(*k)*sizeof(double));
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ ) {
              A[i*LDA+j]=a[j*(*lda)+i];
           }
     }

     if (transb == CblasNoTrans) {
        LDB = *n+1;
        B=(double* )malloc((*k)*LDB*sizeof(double) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ ) {
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }
     else {
        LDB = *k+1;
        B=(double* )malloc(LDB*(*n)*sizeof(double));
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ ) {
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }

     LDC = *n+1;
     C=(double* )malloc((*n)*LDC*sizeof(double));
     for( j=0; j<*n; j++ )
        for( i=0; i<*n; i++ ) {
           C[i*LDC+j]=c[j*(*ldc)+i];
        }
     cblas_dgemmtr( CblasRowMajor, uplo, transa, transb, *n, *k, *alpha, A, LDA,
                  B, LDB, *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*n; i++ ) {
           c[j*(*ldc)+i]=C[i*LDC+j];
        }
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR){
     cblas_dgemmtr( CblasColMajor, uplo, transa, transb, *n, *k, *alpha, a, *lda,
                  b, *ldb, *beta, c, *ldc );
  }
  else
     cblas_dgemmtr( UNDEFINED, uplo, transa, transb, *n, *k, *alpha, a, *lda,
                  b, *ldb, *beta, c, *ldc );
}





void F77_dsymm(CBLAS_INT *layout, char *rtlf, char *uplow, CBLAS_INT *m, CBLAS_INT *n,
              double *alpha, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb,
              double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rtlf_len, FORTRAN_STRLEN uplow_len
#endif
) {

  double *A, *B, *C;
  CBLAS_INT i,j,LDA, LDB, LDC;
  CBLAS_UPLO uplo;
  CBLAS_SIDE side;

  get_uplo_type(uplow,&uplo);
  get_side_type(rtlf,&side);

  if (*layout == TEST_ROW_MJR) {
     if (side == CblasLeft) {
        LDA = *m+1;
        A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );
        for( i=0; i<*m; i++ )
           for( j=0; j<*m; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else{
        LDA = *n+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*n; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     LDB = *n+1;
     B   = ( double* )malloc( (*m)*LDB*sizeof( double ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           B[i*LDB+j]=b[j*(*ldb)+i];
     LDC = *n+1;
     C   = ( double* )malloc( (*m)*LDC*sizeof( double ) );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           C[i*LDC+j]=c[j*(*ldc)+i];
     cblas_dsymm( CblasRowMajor, side, uplo, *m, *n, *alpha, A, LDA, B, LDB,
                  *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dsymm( CblasColMajor, side, uplo, *m, *n, *alpha, a, *lda, b, *ldb,
                  *beta, c, *ldc );
  else
     cblas_dsymm( UNDEFINED, side, uplo, *m, *n, *alpha, a, *lda, b, *ldb,
                  *beta, c, *ldc );
}

void F77_dkymm(CBLAS_INT *layout, char *rtlf, char *uplow, CBLAS_INT *m, CBLAS_INT *n,
              double *alpha, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb,
              double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rtlf_len, FORTRAN_STRLEN uplow_len
#endif
) {

  double *A, *B, *C;
  CBLAS_INT i,j,LDA, LDB, LDC;
  CBLAS_UPLO uplo;
  CBLAS_SIDE side;

  get_uplo_type(uplow,&uplo);
  get_side_type(rtlf,&side);

  if (*layout == TEST_ROW_MJR) {
     if (side == CblasLeft) {
        LDA = *m+1;
        A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );
        for( i=0; i<*m; i++ )
           for( j=0; j<*m; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else{
        LDA = *n+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*n; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     LDB = *n+1;
     B   = ( double* )malloc( (*m)*LDB*sizeof( double ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           B[i*LDB+j]=b[j*(*ldb)+i];
     LDC = *n+1;
     C   = ( double* )malloc( (*m)*LDC*sizeof( double ) );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           C[i*LDC+j]=c[j*(*ldc)+i];
     cblas_dkymm( CblasRowMajor, side, uplo, *m, *n, *alpha, A, LDA, B, LDB,
                  *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dkymm( CblasColMajor, side, uplo, *m, *n, *alpha, a, *lda, b, *ldb,
                  *beta, c, *ldc );
  else
     cblas_dkymm( UNDEFINED, side, uplo, *m, *n, *alpha, a, *lda, b, *ldb,
                  *beta, c, *ldc );
}

void F77_dsyrk(CBLAS_INT *layout, char *uplow, char *transp, CBLAS_INT *n, CBLAS_INT *k,
              double *alpha, double *a, CBLAS_INT *lda,
              double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN uplow_len, FORTRAN_STRLEN transp_len
#endif
) {

  CBLAS_INT i,j,LDA,LDC;
  double *A, *C;
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;

  get_uplo_type(uplow,&uplo);
  get_transpose_type(transp,&trans);

  if (*layout == TEST_ROW_MJR) {
     if (trans == CblasNoTrans) {
        LDA = *k+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else{
        LDA = *n+1;
        A   = ( double* )malloc( (*k)*LDA*sizeof( double ) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     LDC = *n+1;
     C   = ( double* )malloc( (*n)*LDC*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           C[i*LDC+j]=c[j*(*ldc)+i];
     cblas_dsyrk(CblasRowMajor, uplo, trans, *n, *k, *alpha, A, LDA, *beta,
	         C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*n; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dsyrk(CblasColMajor, uplo, trans, *n, *k, *alpha, a, *lda, *beta,
	         c, *ldc );
  else
     cblas_dsyrk(UNDEFINED, uplo, trans, *n, *k, *alpha, a, *lda, *beta,
	         c, *ldc );
}

void F77_dsyr2k(CBLAS_INT *layout, char *uplow, char *transp, CBLAS_INT *n, CBLAS_INT *k,
               double *alpha, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb,
               double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN uplow_len, FORTRAN_STRLEN transp_len
#endif
) {
  CBLAS_INT i,j,LDA,LDB,LDC;
  double *A, *B, *C;
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;

  get_uplo_type(uplow,&uplo);
  get_transpose_type(transp,&trans);

  if (*layout == TEST_ROW_MJR) {
     if (trans == CblasNoTrans) {
        LDA = *k+1;
        LDB = *k+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        B   = ( double* )malloc( (*n)*LDB*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ ) {
              A[i*LDA+j]=a[j*(*lda)+i];
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }
     else {
        LDA = *n+1;
        LDB = *n+1;
        A   = ( double* )malloc( LDA*(*k)*sizeof( double ) );
        B   = ( double* )malloc( LDB*(*k)*sizeof( double ) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ ){
              A[i*LDA+j]=a[j*(*lda)+i];
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }
     LDC = *n+1;
     C   = ( double* )malloc( (*n)*LDC*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           C[i*LDC+j]=c[j*(*ldc)+i];
     cblas_dsyr2k(CblasRowMajor, uplo, trans, *n, *k, *alpha, A, LDA,
		  B, LDB, *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*n; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dsyr2k(CblasColMajor, uplo, trans, *n, *k, *alpha, a, *lda,
		   b, *ldb, *beta, c, *ldc );
  else
     cblas_dsyr2k(UNDEFINED, uplo, trans, *n, *k, *alpha, a, *lda,
		   b, *ldb, *beta, c, *ldc );
}
void F77_dkyr2k(CBLAS_INT *layout, char *uplow, char *transp, CBLAS_INT *n, CBLAS_INT *k,
               double *alpha, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb,
               double *beta, double *c, CBLAS_INT *ldc
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN uplow_len, FORTRAN_STRLEN transp_len
#endif
) {
  CBLAS_INT i,j,LDA,LDB,LDC;
  double *A, *B, *C;
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;

  get_uplo_type(uplow,&uplo);
  get_transpose_type(transp,&trans);

  if (*layout == TEST_ROW_MJR) {
     if (trans == CblasNoTrans) {
        LDA = *k+1;
        LDB = *k+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        B   = ( double* )malloc( (*n)*LDB*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*k; j++ ) {
              A[i*LDA+j]=a[j*(*lda)+i];
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }
     else {
        LDA = *n+1;
        LDB = *n+1;
        A   = ( double* )malloc( LDA*(*k)*sizeof( double ) );
        B   = ( double* )malloc( LDB*(*k)*sizeof( double ) );
        for( i=0; i<*k; i++ )
           for( j=0; j<*n; j++ ){
              A[i*LDA+j]=a[j*(*lda)+i];
              B[i*LDB+j]=b[j*(*ldb)+i];
           }
     }
     LDC = *n+1;
     C   = ( double* )malloc( (*n)*LDC*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           C[i*LDC+j]=c[j*(*ldc)+i];
     cblas_dkyr2k(CblasRowMajor, uplo, trans, *n, *k, *alpha, A, LDA,
		  B, LDB, *beta, C, LDC );
     for( j=0; j<*n; j++ )
        for( i=0; i<*n; i++ )
           c[j*(*ldc)+i]=C[i*LDC+j];
     free(A);
     free(B);
     free(C);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dkyr2k(CblasColMajor, uplo, trans, *n, *k, *alpha, a, *lda,
		   b, *ldb, *beta, c, *ldc );
  else
     cblas_dkyr2k(UNDEFINED, uplo, trans, *n, *k, *alpha, a, *lda,
		   b, *ldb, *beta, c, *ldc );
}
void F77_dtrmm(CBLAS_INT *layout, char *rtlf, char *uplow, char *transp, char *diagn,
              CBLAS_INT *m, CBLAS_INT *n, double *alpha, double *a, CBLAS_INT *lda, double *b,
              CBLAS_INT *ldb
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rtlf_len, FORTRAN_STRLEN uplow_len, FORTRAN_STRLEN transp_len, FORTRAN_STRLEN diag_len
#endif
) {
  CBLAS_INT i,j,LDA,LDB;
  double *A, *B;
  CBLAS_SIDE side;
  CBLAS_DIAG diag;
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;

  get_uplo_type(uplow,&uplo);
  get_transpose_type(transp,&trans);
  get_diag_type(diagn,&diag);
  get_side_type(rtlf,&side);

  if (*layout == TEST_ROW_MJR) {
     if (side == CblasLeft) {
        LDA = *m+1;
        A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );
        for( i=0; i<*m; i++ )
           for( j=0; j<*m; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else{
        LDA = *n+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*n; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     LDB = *n+1;
     B   = ( double* )malloc( (*m)*LDB*sizeof( double ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           B[i*LDB+j]=b[j*(*ldb)+i];
     cblas_dtrmm(CblasRowMajor, side, uplo, trans, diag, *m, *n, *alpha,
		 A, LDA, B, LDB );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           b[j*(*ldb)+i]=B[i*LDB+j];
     free(A);
     free(B);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dtrmm(CblasColMajor, side, uplo, trans, diag, *m, *n, *alpha,
		   a, *lda, b, *ldb);
  else
     cblas_dtrmm(UNDEFINED, side, uplo, trans, diag, *m, *n, *alpha,
		   a, *lda, b, *ldb);
}

void F77_dtrsm(CBLAS_INT *layout, char *rtlf, char *uplow, char *transp, char *diagn,
              CBLAS_INT *m, CBLAS_INT *n, double *alpha, double *a, CBLAS_INT *lda, double *b,
              CBLAS_INT *ldb
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rtlf_len, FORTRAN_STRLEN uplow_len, FORTRAN_STRLEN transp_len, FORTRAN_STRLEN diagn_len
#endif
) {
  CBLAS_INT i,j,LDA,LDB;
  double *A, *B;
  CBLAS_SIDE side;
  CBLAS_DIAG diag;
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;

  get_uplo_type(uplow,&uplo);
  get_transpose_type(transp,&trans);
  get_diag_type(diagn,&diag);
  get_side_type(rtlf,&side);

  if (*layout == TEST_ROW_MJR) {
     if (side == CblasLeft) {
        LDA = *m+1;
        A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );
        for( i=0; i<*m; i++ )
           for( j=0; j<*m; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     else{
        LDA = *n+1;
        A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
        for( i=0; i<*n; i++ )
           for( j=0; j<*n; j++ )
              A[i*LDA+j]=a[j*(*lda)+i];
     }
     LDB = *n+1;
     B   = ( double* )malloc( (*m)*LDB*sizeof( double ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           B[i*LDB+j]=b[j*(*ldb)+i];
     cblas_dtrsm(CblasRowMajor, side, uplo, trans, diag, *m, *n, *alpha,
		 A, LDA, B, LDB );
     for( j=0; j<*n; j++ )
        for( i=0; i<*m; i++ )
           b[j*(*ldb)+i]=B[i*LDB+j];
     free(A);
     free(B);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dtrsm(CblasColMajor, side, uplo, trans, diag, *m, *n, *alpha,
		   a, *lda, b, *ldb);
  else
     cblas_dtrsm(UNDEFINED, side, uplo, trans, diag, *m, *n, *alpha,
		   a, *lda, b, *ldb);
}
