#include <stdio.h>
#include <string.h>
#include "cblas.h"
#include "cblas_test.h"

CBLAS_INT cblas_ok, cblas_lerr, cblas_info;
CBLAS_INT link_xerbla=TRUE;
char *cblas_rout;

#ifdef F77_Char
void F77_xerbla(F77_Char F77_srname, void *vinfo
#else
void F77_xerbla(char *srname, void *vinfo
#endif
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN srname_len
#endif
);

void chkxer(void) {
   extern CBLAS_INT cblas_ok, cblas_lerr, cblas_info;
   extern CBLAS_INT link_xerbla;
   extern char *cblas_rout;
   if (cblas_lerr == 1 ) {
      printf("***** ILLEGAL VALUE OF PARAMETER NUMBER %d NOT DETECTED BY %s *****\n", (int) cblas_info, cblas_rout);
      cblas_ok = 0 ;
   }
   cblas_lerr = 1 ;
}

void F77_d2chke(char *rout
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rout_len
#endif
) {
   char *sf = ( rout ) ;
   double A[2] = {0.0,0.0},
          X[2] = {0.0,0.0},
          Y[2] = {0.0,0.0},
          ALPHA=0.0, BETA=0.0;
   extern CBLAS_INT cblas_info, cblas_lerr, cblas_ok;
   extern char *cblas_rout;

#ifndef HAS_ATTRIBUTE_WEAK_SUPPORT
   #ifdef CBLAS_DLL_IMPORTS
   // Since Windows does not support weak symbols, and the trick below doesn't
   // work for shared libraries on Windows, we skip the xerbla tests here.
   printf("***** WARNING: Skipping xerbla tests since weak symbols are not supported on Windows *****\n");
   return;
   #endif

   if (link_xerbla) /* call these first to link */
   {
      API_SUFFIX(cblas_xerbla)(cblas_info,cblas_rout,"");
      F77_xerbla(cblas_rout,&cblas_info, 1);
   }
#endif

   link_xerbla = 0;
   cblas_ok = TRUE ;
   cblas_lerr = PASSED ;

   if (strncmp( sf,"cblas_dgemv",11)==0) {
      cblas_rout = "cblas_dgemv";
      cblas_info = 1;
      API_SUFFIX(cblas_dgemv)(INVALID_LAYOUT, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, INVALID_TRANSPOSE, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgemv)(CblasColMajor, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();

      cblas_info = 2; RowMajorStrg = TRUE; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, INVALID_TRANSPOSE, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, CblasNoTrans, 0, 2,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgemv)(CblasRowMajor, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dgbmv",11)==0) {
      cblas_rout = "cblas_dgbmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(INVALID_LAYOUT, CblasNoTrans, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, INVALID_TRANSPOSE, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, INVALID, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 0, INVALID, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 0, 0, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 2, 0, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 0, 0, 1, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dgbmv)(CblasColMajor, CblasNoTrans, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, INVALID_TRANSPOSE, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, INVALID, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 0, INVALID, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 0, 0, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 2, 0, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 0, 0, 1, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dgbmv)(CblasRowMajor, CblasNoTrans, 0, 0, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dsymv",11)==0) {
      cblas_rout = "cblas_dsymv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(INVALID_LAYOUT, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(CblasColMajor, INVALID_UPLO, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(CblasColMajor, CblasUpper, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(CblasColMajor, CblasUpper, 2,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsymv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsymv)(CblasRowMajor, INVALID_UPLO, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsymv)(CblasRowMajor, CblasUpper, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsymv)(CblasRowMajor, CblasUpper, 2,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsymv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsymv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dskewsymv",15)==0) {
      cblas_rout = "cblas_dskewsymv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(INVALID_LAYOUT, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(CblasColMajor, INVALID_UPLO, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(CblasColMajor, CblasUpper, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(CblasColMajor, CblasUpper, 2,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsymv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsymv)(CblasRowMajor, INVALID_UPLO, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsymv)(CblasRowMajor, CblasUpper, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsymv)(CblasRowMajor, CblasUpper, 2,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsymv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsymv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dsbmv",11)==0) {
      cblas_rout = "cblas_dsbmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(INVALID_LAYOUT, CblasUpper, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, INVALID_UPLO, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, CblasUpper, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, CblasUpper, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, CblasUpper, 0, 1,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, CblasUpper, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsbmv)(CblasColMajor, CblasUpper, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, INVALID_UPLO, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, CblasUpper, INVALID, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, CblasUpper, 0, INVALID,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, CblasUpper, 0, 1,
                  ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, CblasUpper, 0, 0,
                  ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsbmv)(CblasRowMajor, CblasUpper, 0, 0,
                  ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dspmv",11)==0) {
      cblas_rout = "cblas_dspmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspmv)(INVALID_LAYOUT, CblasUpper, 0,
                  ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspmv)(CblasColMajor, INVALID_UPLO, 0,
                  ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspmv)(CblasColMajor, CblasUpper, INVALID,
                  ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspmv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspmv)(CblasColMajor, CblasUpper, 0,
                  ALPHA, A, X, 1, BETA, Y, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspmv)(CblasRowMajor, INVALID_UPLO, 0,
                  ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspmv)(CblasRowMajor, CblasUpper, INVALID,
                  ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspmv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, X, 0, BETA, Y, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspmv)(CblasRowMajor, CblasUpper, 0,
                  ALPHA, A, X, 1, BETA, Y, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtrmv",11)==0) {
      cblas_rout = "cblas_dtrmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 2, A, 1, X, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 2, A, 1, X, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtbmv",11)==0) {
      cblas_rout = "cblas_dtbmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 1, A, 1, X, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 1, A, 1, X, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtpmv",11)==0) {
      cblas_rout = "cblas_dtpmv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpmv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpmv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpmv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpmv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtrsv",11)==0) {
      cblas_rout = "cblas_dtrsv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 2, A, 1, X, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtrsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 2, A, 1, X, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtrsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, 1, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtbsv",11)==0) {
      cblas_rout = "cblas_dtbsv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 1, A, 1, X, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtbsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, 0, A, 1, X, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, INVALID, A, 1, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 1, A, 1, X, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtbsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, 0, A, 1, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dtpsv",11)==0) {
      cblas_rout = "cblas_dtpsv";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(INVALID_LAYOUT, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(CblasColMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(CblasColMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dtpsv)(CblasColMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 0 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpsv)(CblasRowMajor, INVALID_UPLO, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpsv)(CblasRowMajor, CblasUpper, INVALID_TRANSPOSE,
                  CblasNonUnit, 0, A, X, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  INVALID_DIAG, 0, A, X, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, INVALID, A, X, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dtpsv)(CblasRowMajor, CblasUpper, CblasNoTrans,
                  CblasNonUnit, 0, A, X, 0 );
      chkxer();
   } else if (strncmp( sf,"cblas_dger",10)==0) {
      cblas_rout = "cblas_dger";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(INVALID_LAYOUT, 0, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(CblasColMajor, INVALID, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(CblasColMajor, 0, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(CblasColMajor, 0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(CblasColMajor, 0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dger)(CblasColMajor, 2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dger)(CblasRowMajor, INVALID, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dger)(CblasRowMajor, 0, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dger)(CblasRowMajor, 0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dger)(CblasRowMajor, 0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dger)(CblasRowMajor, 0, 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
   } else if (strncmp( sf,"cblas_dsyr2",11)==0) {
      cblas_rout = "cblas_dsyr2";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(INVALID_LAYOUT, CblasUpper, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr2)(CblasColMajor, CblasUpper, 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr2)(CblasRowMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr2)(CblasRowMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr2)(CblasRowMajor, CblasUpper, 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
   } else if (strncmp( sf,"cblas_dskewsyr2",15)==0) {
      cblas_rout = "cblas_dskewsyr2";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(INVALID_LAYOUT, CblasUpper, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dskewsyr2)(CblasColMajor, CblasUpper, 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsyr2)(CblasRowMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsyr2)(CblasRowMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsyr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsyr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dskewsyr2)(CblasRowMajor, CblasUpper, 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer();
   } else if (strncmp( sf,"cblas_dspr2",11)==0) {
      cblas_rout = "cblas_dspr2";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr2)(INVALID_LAYOUT, CblasUpper, 0, ALPHA, X, 1, Y, 1, A );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr2)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr2)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr2)(CblasColMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspr2)(CblasRowMajor, INVALID_UPLO, 0, ALPHA, X, 1, Y, 1, A );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspr2)(CblasRowMajor, CblasUpper, INVALID, ALPHA, X, 1, Y, 1, A );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 0, Y, 1, A );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dspr2)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 1, Y, 0, A );
      chkxer();
   } else if (strncmp( sf,"cblas_dsyr",10)==0) {
      cblas_rout = "cblas_dsyr";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr)(INVALID_LAYOUT, CblasUpper, 0, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dsyr)(CblasColMajor, CblasUpper, 2, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr)(CblasRowMajor, INVALID_UPLO, 0, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr)(CblasRowMajor, CblasUpper, INVALID, ALPHA, X, 1, A, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr)(CblasRowMajor, CblasUpper, 0, ALPHA, X, 0, A, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      API_SUFFIX(cblas_dsyr)(CblasRowMajor, CblasUpper, 2, ALPHA, X, 1, A, 1 );
      chkxer();
   } else if (strncmp( sf,"cblas_dspr",10)==0) {
      cblas_rout = "cblas_dspr";
      cblas_info = 1; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(INVALID_LAYOUT, CblasUpper, 0, ALPHA, X, 1, A );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, A );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, A );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, A );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, INVALID_UPLO, 0, ALPHA, X, 1, A );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, CblasUpper, INVALID, ALPHA, X, 1, A );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      API_SUFFIX(cblas_dspr)(CblasColMajor, CblasUpper, 0, ALPHA, X, 0, A );
      chkxer();
   }
   if (cblas_ok == TRUE)
       printf(" %-16s PASSED THE TESTS OF ERROR-EXITS\n", cblas_rout);
   else
       printf("******* %s FAILED THE TESTS OF ERROR-EXITS *******\n",cblas_rout);
}
