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

void  F77_c3chke(char *  rout
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rout_len
#endif
) {
   char *sf = ( rout ) ;
   float   A[4]     = {0.0,0.0,0.0,0.0},
           B[4]     = {0.0,0.0,0.0,0.0},
           C[4]     = {0.0,0.0,0.0,0.0},
           ALPHA[2] = {0.0,0.0},
           BETA[2]  = {0.0,0.0},
           RALPHA   = 0.0, RBETA = 0.0;
   extern CBLAS_INT cblas_info, cblas_lerr, cblas_ok;
   extern int RowMajorStrg;
   extern char *cblas_rout;

   cblas_ok = TRUE ;
   cblas_lerr = PASSED ;

#ifndef HAS_ATTRIBUTE_WEAK_SUPPORT
   if (link_xerbla) /* call these first to link */
   {
      cblas_xerbla(cblas_info,cblas_rout,"");
      F77_xerbla(cblas_rout,&cblas_info, 1);
   }
#endif

   link_xerbla = 0;
   if (strncmp( sf,"cblas_cgemmtr"   ,13)==0) {
      cblas_rout = "cblas_cgemmtr"   ;

      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasUpper, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasUpper, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasUpper,CblasTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasUpper, CblasTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasLower, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasLower, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasLower,CblasTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemmtr( INVALID, CblasLower, CblasTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  INVALID, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  INVALID, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, INVALID, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor, CblasUpper,  CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor, CblasUpper, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();


      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();


      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

      /* Row Major */
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 9;  RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_cgemm"   ,11)==0) {
      cblas_rout = "cblas_cgemm"   ;

      cblas_info = 1;
      cblas_cgemm( INVALID,  CblasNoTrans, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemm( INVALID,  CblasNoTrans, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemm( INVALID,  CblasTrans, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_cgemm( INVALID,  CblasTrans, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  INVALID, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  INVALID, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, INVALID, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, INVALID, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_cgemm( CblasColMajor,  CblasTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 9;  RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_cgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
   } else if (strncmp( sf,"cblas_chemm"   ,11)==0) {
            cblas_rout = "cblas_chemm"   ;

      cblas_info = 1;
      cblas_chemm( INVALID,  CblasRight, CblasLower, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  INVALID, CblasUpper, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_chemm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_chemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_csymm"   ,11)==0) {
            cblas_rout = "cblas_csymm"   ;

      cblas_info = 1;
      cblas_csymm( INVALID,  CblasRight, CblasLower, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  INVALID, CblasUpper, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csymm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_ctrmm"   ,11)==0) {
            cblas_rout = "cblas_ctrmm"   ;

      cblas_info = 1;
      cblas_ctrmm( INVALID,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  INVALID, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, INVALID, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, INVALID,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   INVALID, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_ctrsm"   ,11)==0) {
            cblas_rout = "cblas_ctrsm"   ;

      cblas_info = 1;
      cblas_ctrsm( INVALID,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  INVALID, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, INVALID, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, INVALID,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   INVALID, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ctrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ctrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_cherk"   ,11)==0) {
            cblas_rout = "cblas_cherk"   ;

      cblas_info = 1;
      cblas_cherk(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasConjTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasConjTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasConjTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_cherk(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_cherk(CblasColMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_csyrk"   ,11)==0) {
            cblas_rout = "cblas_csyrk"   ;

      cblas_info = 1;
      cblas_csyrk(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_csyrk(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_csyrk(CblasColMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_cher2k"   ,12)==0) {
            cblas_rout = "cblas_cher2k"   ;

      cblas_info = 1;
      cblas_cher2k(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasConjTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasConjTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_cher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_cher2k(CblasColMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_csyr2k"   ,12)==0) {
            cblas_rout = "cblas_csyr2k"   ;

      cblas_info = 1;
      cblas_csyr2k(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_csyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_csyr2k(CblasColMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

   }

   if (cblas_ok == 1 )
       printf(" %-13s PASSED THE TESTS OF ERROR-EXITS\n", cblas_rout);
   else
       printf("***** %s FAILED THE TESTS OF ERROR-EXITS *******\n",cblas_rout);
}
