      SUBROUTINE DERRQP( PATH, NUNIT )
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
*  DERRQP tests the error exits for DGEQPF and DGEQP3.
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
      PARAMETER          ( NMAX = 3 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO, LW
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), TAU( NMAX ), W( 3*NMAX+1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGEQP3, DGEQPF
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
      LW = 3*NMAX + 1
      A( 1, 1 ) = 1.0D+0
      A( 1, 2 ) = 2.0D+0
      A( 2, 2 ) = 3.0D+0
      A( 2, 1 ) = 4.0D+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'QP' ) ) THEN
*
*        Test error exits for QR factorization with pivoting
*
*        DGEQPF
*
         SRNAMT = 'DGEQPF'
         INFOT = 1
         CALL DGEQPF( -1, 0, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'DGEQPF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGEQPF( 0, -1, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'DGEQPF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGEQPF( 2, 0, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'DGEQPF', INFOT, NOUT, LERR, OK )
*
*        DGEQP3
*
         SRNAMT = 'DGEQP3'
         INFOT = 1
         CALL DGEQP3( -1, 0, A, 1, IP, TAU, W, LW, INFO )
         CALL CHKXER( 'DGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGEQP3( 1, -1, A, 1, IP, TAU, W, LW, INFO )
         CALL CHKXER( 'DGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGEQP3( 2, 3, A, 1, IP, TAU, W, LW, INFO )
         CALL CHKXER( 'DGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGEQP3( 2, 2, A, 2, IP, TAU, W, LW-10, INFO )
         CALL CHKXER( 'DGEQP3', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRQP
*
      END
