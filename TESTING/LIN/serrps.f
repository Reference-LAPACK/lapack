      SUBROUTINE SERRPS( PATH, NUNIT )
*
*  -- LAPACK test routine (version 3.1) --
*     Craig Lucas, University of Manchester / NAG Ltd.
*     October, 2008
*
*     .. Scalar Arguments ..
      INTEGER            NUNIT
      CHARACTER*3        PATH
*     ..
*
*  Purpose
*  =======
*
*  SERRPS tests the error exits for the REAL routines
*  for SPSTRF..
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
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      REAL               A( NMAX, NMAX ), WORK( 2*NMAX )
      INTEGER            PIV( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SPSTF2, SPSTRF
*     ..
*     .. Scalars in Common ..
      INTEGER            INFOT, NOUT
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO 110 J = 1, NMAX
         DO 100 I = 1, NMAX
            A( I, J ) = 1.0 / REAL( I+J )
*
  100    CONTINUE
         PIV( J ) = J
         WORK( J ) = 0.
         WORK( NMAX+J ) = 0.
*
  110 CONTINUE
      OK = .TRUE.
*
*
*        Test error exits of the routines that use the Cholesky
*        decomposition of a symmetric positive semidefinite matrix.
*
*        SPSTRF
*
      SRNAMT = 'SPSTRF'
      INFOT = 1
      CALL SPSTRF( '/', 0, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL SPSTRF( 'U', -1, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL SPSTRF( 'U', 2, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTRF', INFOT, NOUT, LERR, OK )
*
*        SPSTF2
*
      SRNAMT = 'SPSTF2'
      INFOT = 1
      CALL SPSTF2( '/', 0, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL SPSTF2( 'U', -1, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL SPSTF2( 'U', 2, A, 1, PIV, 1, -1.0, WORK, INFO )
      CALL CHKXER( 'SPSTF2', INFOT, NOUT, LERR, OK )
*
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRPS
*
      END
