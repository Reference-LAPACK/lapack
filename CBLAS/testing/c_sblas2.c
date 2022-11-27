/*
 *     Written by D.P. Manley, Digital Equipment Corporation.
 *     Prefixed "C_" to BLAS routines and their declarations.
 *
 *     Modified by T. H. Do, 1/23/98, SGI/CRAY Research.
 */
#include <stdlib.h>
#include "cblas.h"
#include "cblas_test.h"

void F77_sgemv(CBLAS_INT *layout, char *transp, CBLAS_INT *m, CBLAS_INT *n, float *alpha,
	       float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx, float *beta,
	       float *y, CBLAS_INT *incy ) {

  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;

  get_transpose_type(transp, &trans);
  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*m)*LDA*sizeof( float ) );
     for( i=0; i<*m; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_sgemv( CblasRowMajor, trans,
		  *m, *n, *alpha, A, LDA, x, *incx, *beta, y, *incy );
     free(A);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_sgemv( CblasColMajor, trans,
		  *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy );
  else
     cblas_sgemv( UNDEFINED, trans,
		  *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy );
}

void F77_sger(CBLAS_INT *layout, CBLAS_INT *m, CBLAS_INT *n, float *alpha, float *x, CBLAS_INT *incx,
	     float *y, CBLAS_INT *incy, float *a, CBLAS_INT *lda ) {

  float *A;
  CBLAS_INT i,j,LDA;

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*m)*LDA*sizeof( float ) );

     for( i=0; i<*m; i++ ) {
       for( j=0; j<*n; j++ )
         A[ LDA*i+j ]=a[ (*lda)*j+i ];
     }

     cblas_sger(CblasRowMajor, *m, *n, *alpha, x, *incx, y, *incy, A, LDA );
     for( i=0; i<*m; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
  }
  else
     cblas_sger( CblasColMajor, *m, *n, *alpha, x, *incx, y, *incy, a, *lda );
}

void F77_strmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx) {
  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*n)*LDA*sizeof( float ) );
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_strmv(CblasRowMajor, uplo, trans, diag, *n, A, LDA, x, *incx);
     free(A);
  }
  else if (*layout == TEST_COL_MJR)
     cblas_strmv(CblasColMajor, uplo, trans, diag, *n, a, *lda, x, *incx);
  else {
     cblas_strmv(UNDEFINED, uplo, trans, diag, *n, a, *lda, x, *incx);
  }
}

void F77_strsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	       CBLAS_INT *n, float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx ) {
  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*n)*LDA*sizeof( float ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_strsv(CblasRowMajor, uplo, trans, diag, *n, A, LDA, x, *incx );
     free(A);
   }
   else
     cblas_strsv(CblasColMajor, uplo, trans, diag, *n, a, *lda, x, *incx );
}
void F77_ssymv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *a,
	      CBLAS_INT *lda, float *x, CBLAS_INT *incx, float *beta, float *y,
	      CBLAS_INT *incy) {
  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*n)*LDA*sizeof( float ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_ssymv(CblasRowMajor, uplo, *n, *alpha, A, LDA, x, *incx,
		 *beta, y, *incy );
     free(A);
   }
   else
     cblas_ssymv(CblasColMajor, uplo, *n, *alpha, a, *lda, x, *incx,
		 *beta, y, *incy );
}

void F77_ssyr(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *x,
	     CBLAS_INT *incx, float *a, CBLAS_INT *lda) {
  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*n)*LDA*sizeof( float ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_ssyr(CblasRowMajor, uplo, *n, *alpha, x, *incx, A, LDA);
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
   }
   else
     cblas_ssyr(CblasColMajor, uplo, *n, *alpha, x, *incx, a, *lda);
}

void F77_ssyr2(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *x,
	     CBLAS_INT *incx, float *y, CBLAS_INT *incy, float *a, CBLAS_INT *lda) {
  float *A;
  CBLAS_INT i,j,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n+1;
     A   = ( float* )malloc( (*n)*LDA*sizeof( float ) );
     for( i=0; i<*n; i++ )
        for( j=0; j<*n; j++ )
           A[ LDA*i+j ]=a[ (*lda)*j+i ];
     cblas_ssyr2(CblasRowMajor, uplo, *n, *alpha, x, *incx, y, *incy, A, LDA);
     for( i=0; i<*n; i++ )
       for( j=0; j<*n; j++ )
         a[ (*lda)*j+i ]=A[ LDA*i+j ];
     free(A);
   }
   else
     cblas_ssyr2(CblasColMajor, uplo, *n, *alpha, x, *incx, y, *incy, a, *lda);
}

void F77_sgbmv(CBLAS_INT *layout, char *transp, CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku,
	       float *alpha, float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx,
	       float *beta, float *y, CBLAS_INT *incy ) {

  float *A;
  CBLAS_INT i,irow,j,jcol,LDA;
  CBLAS_TRANSPOSE trans;

  get_transpose_type(transp, &trans);

  if (*layout == TEST_ROW_MJR) {
     LDA = *ku+*kl+2;
     A   = ( float* )malloc( (*n+*kl)*LDA*sizeof( float ) );
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
     cblas_sgbmv( CblasRowMajor, trans, *m, *n, *kl, *ku, *alpha,
		  A, LDA, x, *incx, *beta, y, *incy );
     free(A);
  }
  else
     cblas_sgbmv( CblasColMajor, trans, *m, *n, *kl, *ku, *alpha,
		  a, *lda, x, *incx, *beta, y, *incy );
}

void F77_stbmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, CBLAS_INT *k, float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx) {
  float *A;
  CBLAS_INT irow, jcol, i, j, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A = ( float* )malloc( (*n+*k)*LDA*sizeof( float ) );
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
     cblas_stbmv(CblasRowMajor, uplo, trans, diag, *n, *k, A, LDA, x, *incx);
     free(A);
   }
   else
     cblas_stbmv(CblasColMajor, uplo, trans, diag, *n, *k, a, *lda, x, *incx);
}

void F77_stbsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, CBLAS_INT *k, float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx) {
  float *A;
  CBLAS_INT irow, jcol, i, j, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A = ( float* )malloc( (*n+*k)*LDA*sizeof( float ) );
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
     cblas_stbsv(CblasRowMajor, uplo, trans, diag, *n, *k, A, LDA, x, *incx);
     free(A);
  }
  else
     cblas_stbsv(CblasColMajor, uplo, trans, diag, *n, *k, a, *lda, x, *incx);
}

void F77_ssbmv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, CBLAS_INT *k, float *alpha,
	      float *a, CBLAS_INT *lda, float *x, CBLAS_INT *incx, float *beta,
	      float *y, CBLAS_INT *incy) {
  float *A;
  CBLAS_INT i,j,irow,jcol,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *k+1;
     A   = ( float* )malloc( (*n+*k)*LDA*sizeof( float ) );
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
     cblas_ssbmv(CblasRowMajor, uplo, *n, *k, *alpha, A, LDA, x, *incx,
		 *beta, y, *incy );
     free(A);
   }
   else
     cblas_ssbmv(CblasColMajor, uplo, *n, *k, *alpha, a, *lda, x, *incx,
		 *beta, y, *incy );
}

void F77_sspmv(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *ap,
	      float *x, CBLAS_INT *incx, float *beta, float *y, CBLAS_INT *incy) {
  float *A,*AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( float* )malloc( LDA*LDA*sizeof( float ) );
     AP  = ( float* )malloc( (((LDA+1)*LDA)/2)*sizeof( float ) );
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
     cblas_sspmv( CblasRowMajor, uplo, *n, *alpha, AP, x, *incx, *beta, y,
		  *incy );
     free(A); free(AP);
  }
  else
     cblas_sspmv( CblasColMajor, uplo, *n, *alpha, ap, x, *incx, *beta, y,
		  *incy );
}

void F77_stpmv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, float *ap, float *x, CBLAS_INT *incx) {
  float *A, *AP;
  CBLAS_INT i, j, k, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( float* )malloc( LDA*LDA*sizeof( float ) );
     AP  = ( float* )malloc( (((LDA+1)*LDA)/2)*sizeof( float ) );
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
     cblas_stpmv( CblasRowMajor, uplo, trans, diag, *n, AP, x, *incx );
     free(A); free(AP);
  }
  else
     cblas_stpmv( CblasColMajor, uplo, trans, diag, *n, ap, x, *incx );
}

void F77_stpsv(CBLAS_INT *layout, char *uplow, char *transp, char *diagn,
	      CBLAS_INT *n, float *ap, float *x, CBLAS_INT *incx) {
  float *A, *AP;
  CBLAS_INT i, j, k, LDA;
  CBLAS_TRANSPOSE trans;
  CBLAS_UPLO uplo;
  CBLAS_DIAG diag;

  get_transpose_type(transp,&trans);
  get_uplo_type(uplow,&uplo);
  get_diag_type(diagn,&diag);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( float* )malloc( LDA*LDA*sizeof( float ) );
     AP  = ( float* )malloc( (((LDA+1)*LDA)/2)*sizeof( float ) );
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
     cblas_stpsv( CblasRowMajor, uplo, trans, diag, *n, AP, x, *incx );
     free(A); free(AP);
  }
  else
     cblas_stpsv( CblasColMajor, uplo, trans, diag, *n, ap, x, *incx );
}

void F77_sspr(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *x,
	     CBLAS_INT *incx, float *ap ){
  float *A, *AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( float* )malloc( LDA*LDA*sizeof( float ) );
     AP  = ( float* )malloc( (((LDA+1)*LDA)/2)*sizeof( float ) );
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
     cblas_sspr( CblasRowMajor, uplo, *n, *alpha, x, *incx, AP );
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
     free(A); free(AP);
  }
  else
     cblas_sspr( CblasColMajor, uplo, *n, *alpha, x, *incx, ap );
}

void F77_sspr2(CBLAS_INT *layout, char *uplow, CBLAS_INT *n, float *alpha, float *x,
	     CBLAS_INT *incx, float *y, CBLAS_INT *incy, float *ap ){
  float *A, *AP;
  CBLAS_INT i,j,k,LDA;
  CBLAS_UPLO uplo;

  get_uplo_type(uplow,&uplo);

  if (*layout == TEST_ROW_MJR) {
     LDA = *n;
     A   = ( float* )malloc( LDA*LDA*sizeof( float ) );
     AP  = ( float* )malloc( (((LDA+1)*LDA)/2)*sizeof( float ) );
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
     cblas_sspr2( CblasRowMajor, uplo, *n, *alpha, x, *incx, y, *incy, AP );
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
     cblas_sspr2( CblasColMajor, uplo, *n, *alpha, x, *incx, y, *incy, ap );
}
