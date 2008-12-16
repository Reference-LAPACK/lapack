      SUBROUTINE SERRLS( PATH, NUNIT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  SERRLS tests the error exits for the REAL least squares
*  driver routines (SGELS, SGELSS, SGELSX, SGELSY, SGELSD).
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name for the routines to be tested.
*
*  NUNIT   (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO, IRNK
      REAL               RCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      REAL               A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ),
     $                   W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SGELS, SGELSD, SGELSS, SGELSX,
     $                   SGELSY
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
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.0E+0
      A( 1, 2 ) = 2.0E+0
      A( 2, 2 ) = 3.0E+0
      A( 2, 1 ) = 4.0E+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        Test error exits for the least squares driver routines.
*
*        SGELS
*
         SRNAMT = 'SGELS '
         INFOT = 1
         CALL SGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'SGELS ', INFOT, NOUT, LERR, OK )
*
*        SGELSS
*
         SRNAMT = 'SGELSS'
         INFOT = 1
         CALL SGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'SGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'SGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'SGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'SGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'SGELSS', INFOT, NOUT, LERR, OK )
*
*        SGELSX
*
         SRNAMT = 'SGELSX'
         INFOT = 1
         CALL SGELSX( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'SGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGELSX( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'SGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGELSX( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'SGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SGELSX( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'SGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SGELSX( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'SGELSX', INFOT, NOUT, LERR, OK )
*
*        SGELSY
*
         SRNAMT = 'SGELSY'
         INFOT = 1
         CALL SGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL SGELSY( 2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'SGELSY', INFOT, NOUT, LERR, OK )
*
*        SGELSD
*
         SRNAMT = 'SGELSD'
         INFOT = 1
         CALL SGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10,
     $                IP, INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10,
     $                IP, INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10,
     $                IP, INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10,
     $                IP, INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10,
     $                IP, INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL SGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP,
     $                INFO )
         CALL CHKXER( 'SGELSD', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRLS
*
      END
