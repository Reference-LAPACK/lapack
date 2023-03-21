/*
 *     Written by D.P. Manley, Digital Equipment Corporation.
 *     Prefixed "C_" to BLAS routines and their declarations.
 *
 *     Modified by T. H. Do, 1/23/98, SGI/CRAY Research.
 */
#include <stdlib.h>
#include "cblas.h"
#include "cblas_test.h"

void F77_dgemv(CBLAS_INT *layout, char *transp, CBLAS_INT *m, CBLAS_INT *n, double *alpha,
	       double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx, double *beta,
	       double *y, CBLAS_INT *incy ) {

  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;

  get_transpose_type(transp, &trans);
  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dgemv( CblasRowMajor, trans,
		  *m, *n, *alpha, A, LDA, x, *incx, *beta, y, *incy );
     free(A);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dgemv( CblasColMajor, trans,
		  *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy );
  else
     cblas_dgemv( UNDEFINED, trans,
		  *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy );
}

void F77_dger(CBLAS_INT *layout, CBLAS_INT *m, CBLAS_INT *n, double *alpha, double *x, CBLAS_INT *incx,
	     double *y, CBLAS_INT *incy, double *a, CBLAS_INT *lda ) {

  double *A;
  CBLAS_INT i,j,LDA;

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*m)*LDA*sizeof( double ) );

     for( i=0; i<*m; i++ ) {
       for( j=0; j<*n; j++ )
         A[ LDA*i+j ]=a[ (*lda)*j+i ];
     }

     cblas_dger(CblasRowMajor, *m, *n, *alpha, x, *incx, y, *incy, A, LDA );
     for( i=0; i<*m; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
  }
  else
     cblas_dger( CblasColMajor, *m, *n, *alpha, x, *incx, y, *incy, a, *lda );
}

void F77_dtrmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx) {
  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dtrmv(CblasRowMajor, uplo, trans, diag, *n, A, LDA, x, *incx);
     free(A);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_dtrmv(CblasColMajor, uplo, trans, diag, *n, a, *lda, x, *incx);
  else {
     cblas_dtrmv(UNDEFINED, uplo, trans, diag, *n, a, *lda, x, *incx);
  }
}

void F77_dtrsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	       CBLAS_INT *n, double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx ) {
  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dtrsv(CblasRowMajor, uplo, trans, diag, *n, A, LDA, x, *incx );
     free(A);
   }
   else
     cblas_dtrsv(CblasColMajor, uplo, trans, diag, *n, a, *lda, x, *incx );
}
void F77_dsymv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *a,
	      CBLAS_INT *lda, double *x, CBLAS_INT *incx, double *beta, double *y,
	      CBLAS_INT *incy) {
  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dsymv(CblasRowMajor, uplo, *n, *alpha, A, LDA, x, *incx,
		 *beta, y, *incy );
     free(A);
   }
   else
     cblas_dsymv(CblasColMajor, uplo, *n, *alpha, a, *lda, x, *incx,
		 *beta, y, *incy );
}

void F77_dsyr(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *x,
	     CBLAS_INT *incx, double *a, CBLAS_INT *lda) {
  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dsyr(CblasRowMajor, uplo, *n, *alpha, x, *incx, A, LDA);
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
   }
   else
     cblas_dsyr(CblasColMajor, uplo, *n, *alpha, x, *incx, a, *lda);
}

void F77_dsyr2(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *x,
	     CBLAS_INT *incx, double *y, CBLAS_INT *incy, double *a, CBLAS_INT *lda) {
  double *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( double* )malloc( (*n)*LDA*sizeof( double ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_dsyr2(CblasRowMajor, uplo, *n, *alpha, x, *incx, y, *incy, A, LDA);
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
   }
   else
     cblas_dsyr2(CblasColMajor, uplo, *n, *alpha, x, *incx, y, *incy, a, *lda);
}

void F77_dgbmv(CBLAS_INT *layout, char *transp, CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku,
	       double *alpha, double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx,
	       double *beta, double *y, CBLAS_INT *incy ) {

  double *A;
  CBLAS_INT i,irow,j,jcol,LDA;
  CBLAS_TRANSPOSE trans;

  get_transpose_type(transp, &trans);

  if (*layout == TEST_ROW_MJR) {
     LDA = *ku+*kl+2;
     A   = ( double* )malloc( (*n+*kl)*LDA*sizeof( double ) );
     for( i=0; i<*ku; i++ ){
        irow=*ku+*kl-i;
        jcol=(*ku)-i;
        for( j=jcol; j<*n; j++ )
           A[ LDA*(j-jcol)+irow ]=a[ (*lda)*j+i ];
     }
     i=*ku;
     irow=*ku+*kl-i;
     for( j=0; j<*n; j++ )
        A[ LDA*j+irow ]=a[ (*lda)*j+i ];
     for( i=*ku+1; i<*ku+*kl+1; i++ ){
        irow=*ku+*kl-i;
        jcol=i-(*ku);
        for( j=jcol; j<(*n+*kl); j++ )
           A[ LDA*j+irow ]=a[ (*lda)*(j-jcol)+i ];
     }
     cblas_dgbmv( CblasRowMajor, trans, *m, *n, *kl, *ku, *alpha,
		  A, LDA, x, *incx, *beta, y, *incy );
     free(A);
  }
  else
     cblas_dgbmv( CblasColMajor, trans, *m, *n, *kl, *ku, *alpha,
		  a, *lda, x, *incx, *beta, y, *incy );
}

void F77_dtbmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, CBLAS_INT *k, double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx) {
  double *A;
  CBLAS_INT irow, jcol, i, j, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A = ( double* )malloc( (*n+*k)*LDA*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( i=0; i<*k; i++ ){
           irow=*k-i;
           jcol=(*k)-i;
           for( j=jcol; j<*n; j++ )
              A[ LDA*(j-jcol)+irow ]=a[ (*lda)*j+i ];
        }
        i=*k;
        irow=*k-i;
        for( j=0; j<*n; j++ )
           A[ LDA*j+irow ]=a[ (*lda)*j+i ];
     }
     else {
       i=0;
       irow=*k-i;
       for( j=0; j<*n; j++ )
          A[ LDA*j+irow ]=a[ (*lda)*j+i ];
       for( i=1; i<*k+1; i++ ){
          irow=*k-i;
          jcol=i;
          for( j=jcol; j<(*n+*k); j++ )
             A[ LDA*j+irow ]=a[ (*lda)*(j-jcol)+i ];
       }
     }
     cblas_dtbmv(CblasRowMajor, uplo, trans, diag, *n, *k, A, LDA, x, *incx);
     free(A);
   }
   else
     cblas_dtbmv(CblasColMajor, uplo, trans, diag, *n, *k, a, *lda, x, *incx);
}

void F77_dtbsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, CBLAS_INT *k, double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx) {
  double *A;
  CBLAS_INT irow, jcol, i, j, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A = ( double* )malloc( (*n+*k)*LDA*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( i=0; i<*k; i++ ){
        irow=*k-i;
        jcol=(*k)-i;
        for( j=jcol; j<*n; j++ )
           A[ LDA*(j-jcol)+irow ]=a[ (*lda)*j+i ];
        }
        i=*k;
        irow=*k-i;
        for( j=0; j<*n; j++ )
           A[ LDA*j+irow ]=a[ (*lda)*j+i ];
     }
     else {
        i=0;
        irow=*k-i;
        for( j=0; j<*n; j++ )
           A[ LDA*j+irow ]=a[ (*lda)*j+i ];
        for( i=1; i<*k+1; i++ ){
           irow=*k-i;
           jcol=i;
           for( j=jcol; j<(*n+*k); j++ )
              A[ LDA*j+irow ]=a[ (*lda)*(j-jcol)+i ];
        }
     }
     cblas_dtbsv(CblasRowMajor, uplo, trans, diag, *n, *k, A, LDA, x, *incx);
     free(A);
  }
  else
     cblas_dtbsv(CblasColMajor, uplo, trans, diag, *n, *k, a, *lda, x, *incx);
}

void F77_dsbmv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, CBLAS_INT *k, double *alpha,
	      double *a, CBLAS_INT *lda, double *x, CBLAS_INT *incx, double *beta,
	      double *y, CBLAS_INT *incy) {
  double *A;
  CBLAS_INT i,j,irow,jcol,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A   = ( double* )malloc( (*n+*k)*LDA*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( i=0; i<*k; i++ ){
           irow=*k-i;
           jcol=(*k)-i;
           for( j=jcol; j<*n; j++ )
        A[ LDA*(j-jcol)+irow ]=a[ (*lda)*j+i ];
        }
        i=*k;
        irow=*k-i;
        for( j=0; j<*n; j++ )
           A[ LDA*j+irow ]=a[ (*lda)*j+i ];
     }
     else {
        i=0;
        irow=*k-i;
        for( j=0; j<*n; j++ )
           A[ LDA*j+irow ]=a[ (*lda)*j+i ];
        for( i=1; i<*k+1; i++ ){
           irow=*k-i;
           jcol=i;
           for( j=jcol; j<(*n+*k); j++ )
              A[ LDA*j+irow ]=a[ (*lda)*(j-jcol)+i ];
        }
     }
     cblas_dsbmv(CblasRowMajor, uplo, *n, *k, *alpha, A, LDA, x, *incx,
		 *beta, y, *incy );
     free(A);
   }
   else
     cblas_dsbmv(CblasColMajor, uplo, *n, *k, *alpha, a, *lda, x, *incx,
		 *beta, y, *incy );
}

void F77_dspmv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *ap,
	      double *x, CBLAS_INT *incx, double *beta, double *y, CBLAS_INT *incy) {
  double *A,*AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( double* )malloc( LDA*LDA*sizeof( double ) );
     AP  = ( double* )malloc( (((LDA+1)*LDA)/2)*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     else {
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     cblas_dspmv( CblasRowMajor, uplo, *n, *alpha, AP, x, *incx, *beta, y,
		  *incy );
     free(A);
     free(AP);
  }
  else
     cblas_dspmv( CblasColMajor, uplo, *n, *alpha, ap, x, *incx, *beta, y,
		  *incy );
}

void F77_dtpmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, double *ap, double *x, CBLAS_INT *incx) {
  double *A, *AP;
  CBLAS_INT i, j, k, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( double* )malloc( LDA*LDA*sizeof( double ) );
     AP  = ( double* )malloc( (((LDA+1)*LDA)/2)*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     else {
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     cblas_dtpmv( CblasRowMajor, uplo, trans, diag, *n, AP, x, *incx );
     free(A);
     free(AP);
  }
  else
     cblas_dtpmv( CblasColMajor, uplo, trans, diag, *n, ap, x, *incx );
}

void F77_dtpsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, double *ap, double *x, CBLAS_INT *incx) {
  double *A, *AP;
  CBLAS_INT i, j, k, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( double* )malloc( LDA*LDA*sizeof( double ) );
     AP  = ( double* )malloc( (((LDA+1)*LDA)/2)*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];

     }
     else {
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     cblas_dtpsv( CblasRowMajor, uplo, trans, diag, *n, AP, x, *incx );
     free(A);
     free(AP);
  }
  else
     cblas_dtpsv( CblasColMajor, uplo, trans, diag, *n, ap, x, *incx );
}

void F77_dspr(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *x,
	     CBLAS_INT *incx, double *ap ){
  double *A, *AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( double* )malloc( LDA*LDA*sizeof( double ) );
     AP  = ( double* )malloc( (((LDA+1)*LDA)/2)*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     else {
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     cblas_dspr( CblasRowMajor, uplo, *n, *alpha, x, *incx, AP );
     if (uplo == CblasUpper) {
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              A[ LDA*i+j ]=AP[ k ];
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              ap[ k ]=A[ LDA*i+j ];
     }
     else {
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              A[ LDA*i+j ]=AP[ k ];
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              ap[ k ]=A[ LDA*i+j ];
     }
     free(A);
     free(AP);
  }
  else
     cblas_dspr( CblasColMajor, uplo, *n, *alpha, x, *incx, ap );
}

void F77_dspr2(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, double *alpha, double *x,
	     CBLAS_INT *incx, double *y, CBLAS_INT *incy, double *ap ){
  double *A, *AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( double* )malloc( LDA*LDA*sizeof( double ) );
     AP  = ( double* )malloc( (((LDA+1)*LDA)/2)*sizeof( double ) );
     if (uplo == CblasUpper) {
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     else {
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              A[ LDA*i+j ]=ap[ k ];
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              AP[ k ]=A[ LDA*i+j ];
     }
     cblas_dspr2( CblasRowMajor, uplo, *n, *alpha, x, *incx, y, *incy, AP );
     if (uplo == CblasUpper) {
        for( i=0, k=0; i<*n; i++ )
           for( j=i; j<*n; j++, k++ )
              A[ LDA*i+j ]=AP[ k ];
        for( j=0, k=0; j<*n; j++ )
           for( i=0; i<j+1; i++, k++ )
              ap[ k ]=A[ LDA*i+j ];
     }
     else {
        for( i=0, k=0; i<*n; i++ )
           for( j=0; j<i+1; j++, k++ )
              A[ LDA*i+j ]=AP[ k ];
        for( j=0, k=0; j<*n; j++ )
           for( i=j; i<*n; i++, k++ )
              ap[ k ]=A[ LDA*i+j ];
     }
     free(A);
     free(AP);
  }
  else
     cblas_dspr2( CblasColMajor, uplo, *n, *alpha, x, *incx, y, *incy, ap );
}
