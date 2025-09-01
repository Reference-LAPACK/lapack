*> \brief \b CERREC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CERREC( PATH, NUNIT )
*
*       .. Scalar Arguments ..
*       CHARACTER*3        PATH
*       INTEGER            NUNIT
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CERREC tests the error exits for the routines for eigen- condition
*> estimation for REAL matrices:
*>    CTRSYL, CTRSYL3, CTREXC, CTRSNA and CTRSEN.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name for the routines to be tested.
*> \endverbatim
*>
*> \param[in] NUNIT
*> \verbatim
*>          NUNIT is INTEGER
*>          The unit number for output.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CERREC( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LW
      PARAMETER          ( NMAX = 4, LW = NMAX*( NMAX+2 ) )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IFST, ILST, INFO, J, M, NT
      REAL               SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            SEL( NMAX )
      REAL               RW( LW ), S( NMAX ), SEP( NMAX ), SWORK( NMAX )
      COMPLEX            A( NMAX, NMAX ), B( NMAX, NMAX ),
     $                   C( NMAX, NMAX ), WORK( LW ), X( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, CTREXC, CTRSEN, CTRSNA, CTRSYL, CTRSYL3
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      OK = .TRUE.
      NT = 0
*
*     Initialize A, B and SEL
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = ZERO
            B( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, NMAX
         A( I, I ) = ONE
         SEL( I ) = .TRUE.
   30 CONTINUE
*
*     Test CTRSYL
*
      SRNAMT = 'CTRSYL'
      INFOT = 1
      CALL CTRSYL( 'X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTRSYL( 'N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTRSYL( 'N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTRSYL( 'N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTRSYL( 'N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTRSYL( 'N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTRSYL( 'N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL CTRSYL( 'N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'CTRSYL', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Test CTRSYL3
*
      SRNAMT = 'CTRSYL3'
      INFOT = 1
      CALL CTRSYL3( 'X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTRSYL3( 'N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTRSYL3( 'N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTRSYL3( 'N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTRSYL3( 'N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTRSYL3( 'N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTRSYL3( 'N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL CTRSYL3( 'N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE,
     $              SWORK, NMAX, INFO )
      CALL CHKXER( 'CTRSYL3', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Test CTREXC
*
      SRNAMT = 'CTREXC'
      IFST = 1
      ILST = 1
      INFOT = 1
      CALL CTREXC( 'X', 1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTREXC( 'N', -1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 4
      ILST = 2
      CALL CTREXC( 'N', 2, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL CTREXC( 'V', 2, A, 2, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 7
      IFST = 0
      ILST = 1
      CALL CTREXC( 'V', 1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 7
      IFST = 2
      CALL CTREXC( 'V', 1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 8
      IFST = 1
      ILST = 0
      CALL CTREXC( 'V', 1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      INFOT = 8
      ILST = 2
      CALL CTREXC( 'V', 1, A, 1, B, 1, IFST, ILST, INFO )
      CALL CHKXER( 'CTREXC', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Test CTRSNA
*
      SRNAMT = 'CTRSNA'
      INFOT = 1
      CALL CTRSNA( 'X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTRSNA( 'B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTRSNA( 'B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL CTRSNA( 'V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M,
     $             WORK, 2, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL CTRSNA( 'B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M,
     $             WORK, 2, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL CTRSNA( 'B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M,
     $             WORK, 2, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL CTRSNA( 'B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL CTRSNA( 'B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 16
      CALL CTRSNA( 'B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M,
     $             WORK, 1, RW, INFO )
      CALL CHKXER( 'CTRSNA', INFOT, NOUT, LERR, OK )
      NT = NT + 9
*
*     Test CTRSEN
*
      SEL( 1 ) = .FALSE.
      SRNAMT = 'CTRSEN'
      INFOT = 1
      CALL CTRSEN( 'X', 'N', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 1, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTRSEN( 'N', 'X', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 1, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTRSEN( 'N', 'N', SEL, -1, A, 1, B, 1, X, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL CTRSEN( 'N', 'N', SEL, 2, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 2, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL CTRSEN( 'N', 'V', SEL, 2, A, 2, B, 1, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 1, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 14
      CALL CTRSEN( 'N', 'V', SEL, 2, A, 2, B, 2, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 0, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 14
      CALL CTRSEN( 'E', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 1, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 14
      CALL CTRSEN( 'V', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ),
     $             WORK, 3, INFO )
      CALL CHKXER( 'CTRSEN', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (',
     $      I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
      RETURN
*
*     End of CERREC
*
      END
