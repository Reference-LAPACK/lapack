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

void F77_z3chke(char *rout
#ifdef BLAS_FORTRAN_STRLEN_END
  , FORTRAN_STRLEN rout_len
#endif
) {
   char *sf = ( rout ) ;
   double  A[4]     = {0.0,0.0,0.0,0.0},
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
   if (strncmp( sf,"cblas_zgemmtr"   ,13)==0) {
      cblas_rout = "cblas_zgemmtr"   ;

      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasUpper, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasUpper, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasUpper,CblasTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasUpper, CblasTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasLower, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasLower, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasLower,CblasTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemmtr( INVALID, CblasLower, CblasTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  INVALID, CblasNoTrans, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  INVALID, CblasNoTrans, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, INVALID, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor, CblasUpper,  CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor, CblasUpper, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();


      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();


      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

      /* Row Major */
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 9;  RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasColMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasNoTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemmtr( CblasRowMajor,  CblasUpper, CblasTrans, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zgemm"   ,11)==0) {
      cblas_rout = "cblas_zgemm"   ;

      cblas_info = 1;
      cblas_zgemm( INVALID,  CblasNoTrans, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemm( INVALID,  CblasNoTrans, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemm( INVALID,  CblasTrans, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 1;
      cblas_zgemm( INVALID,  CblasTrans, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  INVALID, CblasNoTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  INVALID, CblasTrans, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, INVALID, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, INVALID, 0, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasNoTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = FALSE;
      cblas_zgemm( CblasColMajor,  CblasTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 9;  RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 9; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, 2, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasNoTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 14; RowMajorStrg = TRUE;
      cblas_zgemm( CblasRowMajor,  CblasTrans, CblasTrans, 0, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zhemm"   ,11)==0) {
            cblas_rout = "cblas_zhemm"   ;

      cblas_info = 1;
      cblas_zhemm( INVALID,  CblasRight, CblasLower, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  INVALID, CblasUpper, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zhemm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zhemm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zsymm"   ,11)==0) {
            cblas_rout = "cblas_zsymm"   ;

      cblas_info = 1;
      cblas_zsymm( INVALID,  CblasRight, CblasLower, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  INVALID, CblasUpper, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, INVALID, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsymm( CblasColMajor,  CblasRight, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasUpper, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasLower, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasUpper, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasLower, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasUpper, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasLower, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasUpper, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasUpper, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasLeft, CblasLower, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsymm( CblasRowMajor,  CblasRight, CblasLower, 0, 2,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_ztrmm"   ,11)==0) {
            cblas_rout = "cblas_ztrmm"   ;

      cblas_info = 1;
      cblas_ztrmm( INVALID,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  INVALID, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, INVALID, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, INVALID,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   INVALID, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrmm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrmm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_ztrsm"   ,11)==0) {
            cblas_rout = "cblas_ztrsm"   ;

      cblas_info = 1;
      cblas_ztrsm( INVALID,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  INVALID, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, INVALID, CblasNoTrans,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, INVALID,
                   CblasNonUnit, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   INVALID, 0, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = FALSE;
      cblas_ztrsm( CblasColMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 6; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, INVALID, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 7; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, INVALID, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 2, 0, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 2 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasUpper, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasLeft, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 1, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasNoTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();
      cblas_info = 12; RowMajorStrg = TRUE;
      cblas_ztrsm( CblasRowMajor,  CblasRight, CblasLower, CblasTrans,
                   CblasNonUnit, 0, 2, ALPHA, A, 2, B, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zherk"   ,11)==0) {
            cblas_rout = "cblas_zherk"   ;

      cblas_info = 1;
      cblas_zherk(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasTrans, 0, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasConjTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasConjTrans, INVALID, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasConjTrans, 0, INVALID,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zherk(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasUpper, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  RALPHA, A, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zherk(CblasColMajor,  CblasLower, CblasConjTrans, 2, 0,
                  RALPHA, A, 1, RBETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zsyrk"   ,11)==0) {
            cblas_rout = "cblas_zsyrk"   ;

      cblas_info = 1;
      cblas_zsyrk(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasTrans, INVALID, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasTrans, 0, INVALID,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = TRUE;
      cblas_zsyrk(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasUpper, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                  ALPHA, A, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 11; RowMajorStrg = FALSE;
      cblas_zsyrk(CblasColMajor,  CblasLower, CblasTrans, 2, 0,
                  ALPHA, A, 1, BETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zher2k"   ,12)==0) {
            cblas_rout = "cblas_zher2k"   ;

      cblas_info = 1;
      cblas_zher2k(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasConjTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasConjTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasConjTrans, 0, 2,
                   ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zher2k(CblasRowMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasUpper, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zher2k(CblasColMajor,  CblasLower, CblasConjTrans, 2, 0,
                   ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer();

   } else if (strncmp( sf,"cblas_zsyr2k"   ,12)==0) {
            cblas_rout = "cblas_zsyr2k"   ;

      cblas_info = 1;
      cblas_zsyr2k(INVALID,  CblasUpper, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 2; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  INVALID, CblasNoTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 3; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasConjTrans, 0, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasNoTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 4; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasTrans, INVALID, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 5; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, INVALID,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 8; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer();
      cblas_info = 10; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasTrans, 0, 2,
                   ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = TRUE;
      cblas_zsyr2k(CblasRowMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasUpper, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasNoTrans, 2, 0,
                   ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer();
      cblas_info = 13; RowMajorStrg = FALSE;
      cblas_zsyr2k(CblasColMajor,  CblasLower, CblasTrans, 2, 0,
                   ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer();

   }

   if (cblas_ok == 1 )
       printf(" %-13s PASSED THE TESTS OF ERROR-EXITS\n", cblas_rout);
   else
       printf("***** %s FAILED THE TESTS OF ERROR-EXITS *******\n",cblas_rout);
}
