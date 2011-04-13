      SUBROUTINE DERRTZ( PATH, NUNIT )
*
*  -- LAPACK test routine (version 3.3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  DERRTZ tests the error exits for DTZRQF and STZRZF.
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
      INTEGER            INFO
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   A( NMAX, NMAX ), TAU( NMAX ), W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DTZRQF, DTZRZF
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
      A( 1, 1 ) = 1.D+0
      A( 1, 2 ) = 2.D+0
      A( 2, 2 ) = 3.D+0
      A( 2, 1 ) = 4.D+0
      W( 1 ) = 0.0D+0
      W( 2 ) = 0.0D+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'TZ' ) ) THEN
*
*        Test error exits for the trapezoidal routines.
*
*        DTZRQF
*
         SRNAMT = 'DTZRQF'
         INFOT = 1
         CALL DTZRQF( -1, 0, A, 1, TAU, INFO )
         CALL CHKXER( 'DTZRQF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTZRQF( 1, 0, A, 1, TAU, INFO )
         CALL CHKXER( 'DTZRQF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTZRQF( 2, 2, A, 1, TAU, INFO )
         CALL CHKXER( 'DTZRQF', INFOT, NOUT, LERR, OK )
*
*        DTZRZF
*
         SRNAMT = 'DTZRZF'
         INFOT = 1
         CALL DTZRZF( -1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTZRZF( 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTZRZF( 2, 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTZRZF( 2, 2, A, 2, TAU, W, 0, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTZRZF( 2, 3, A, 2, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRTZ
*
      END
