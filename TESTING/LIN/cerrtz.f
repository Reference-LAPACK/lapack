      SUBROUTINE CERRTZ( PATH, NUNIT )
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
*  CERRTZ tests the error exits for CTZRQF and CTZRZF.
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
      COMPLEX            A( NMAX, NMAX ), TAU( NMAX ), W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, CTZRQF, CTZRZF
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
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = CMPLX( 1.E+0, -1.E+0 )
      A( 1, 2 ) = CMPLX( 2.E+0, -2.E+0 )
      A( 2, 2 ) = CMPLX( 3.E+0, -3.E+0 )
      A( 2, 1 ) = CMPLX( 4.E+0, -4.E+0 )
      W( 1 ) = CMPLX( 0.E+0, 0.E+0 )
      W( 2 ) = CMPLX( 0.E+0, 0.E+0 )
      OK = .TRUE.
*
*     Test error exits for the trapezoidal routines.
*
      WRITE( NOUT, FMT = * )
      IF( LSAMEN( 2, C2, 'TZ' ) ) THEN
*
*        CTZRQF
*
         SRNAMT = 'CTZRQF'
         INFOT = 1
         CALL CTZRQF( -1, 0, A, 1, TAU, INFO )
         CALL CHKXER( 'CTZRQF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTZRQF( 1, 0, A, 1, TAU, INFO )
         CALL CHKXER( 'CTZRQF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTZRQF( 2, 2, A, 1, TAU, INFO )
         CALL CHKXER( 'CTZRQF', INFOT, NOUT, LERR, OK )
*
*        CTZRZF
*
         SRNAMT = 'CTZRZF'
         INFOT = 1
         CALL CTZRZF( -1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTZRZF( 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTZRZF( 2, 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CTZRZF( 2, 2, A, 2, TAU, W, 0, INFO )
         CALL CHKXER( 'CTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CTZRZF( 2, 3, A, 2, TAU, W, 1, INFO )
         CALL CHKXER( 'CTZRZF', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRTZ
*
      END
