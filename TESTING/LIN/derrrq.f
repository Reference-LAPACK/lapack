      SUBROUTINE DERRRQ( PATH, NUNIT )
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
*  DERRRQ tests the error exits for the DOUBLE PRECISION routines
*  that use the RQ decomposition of a general matrix.
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
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   W( NMAX ), X( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGERQ2, DGERQF, DGERQS, DORGR2,
     $                   DORGRQ, DORMR2, DORMRQ
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
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      OK = .TRUE.
*
*     Error exits for RQ factorization
*
*     DGERQF
*
      SRNAMT = 'DGERQF'
      INFOT = 1
      CALL DGERQF( -1, 0, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGERQF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGERQF( 0, -1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGERQF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGERQF( 2, 1, A, 1, B, W, 2, INFO )
      CALL CHKXER( 'DGERQF', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DGERQF( 2, 1, A, 2, B, W, 1, INFO )
      CALL CHKXER( 'DGERQF', INFOT, NOUT, LERR, OK )
*
*     DGERQ2
*
      SRNAMT = 'DGERQ2'
      INFOT = 1
      CALL DGERQ2( -1, 0, A, 1, B, W, INFO )
      CALL CHKXER( 'DGERQ2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGERQ2( 0, -1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGERQ2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGERQ2( 2, 1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGERQ2', INFOT, NOUT, LERR, OK )
*
*     DGERQS
*
      SRNAMT = 'DGERQS'
      INFOT = 1
      CALL DGERQS( -1, 0, 0, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGERQS( 0, -1, 0, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGERQS( 2, 1, 0, A, 2, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DGERQS( 0, 0, -1, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DGERQS( 2, 2, 0, A, 1, X, B, 2, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DGERQS( 2, 2, 0, A, 2, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DGERQS( 1, 1, 2, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGERQS', INFOT, NOUT, LERR, OK )
*
*     DORGRQ
*
      SRNAMT = 'DORGRQ'
      INFOT = 1
      CALL DORGRQ( -1, 0, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGRQ( 0, -1, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGRQ( 2, 1, 0, A, 2, X, W, 2, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGRQ( 0, 0, -1, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGRQ( 1, 2, 2, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORGRQ( 2, 2, 0, A, 1, X, W, 2, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DORGRQ( 2, 2, 0, A, 2, X, W, 1, INFO )
      CALL CHKXER( 'DORGRQ', INFOT, NOUT, LERR, OK )
*
*     DORGR2
*
      SRNAMT = 'DORGR2'
      INFOT = 1
      CALL DORGR2( -1, 0, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGR2( 0, -1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGR2( 2, 1, 0, A, 2, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGR2( 0, 0, -1, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGR2( 1, 2, 2, A, 2, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORGR2( 2, 2, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGR2', INFOT, NOUT, LERR, OK )
*
*     DORMRQ
*
      SRNAMT = 'DORMRQ'
      INFOT = 1
      CALL DORMRQ( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORMRQ( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORMRQ( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORMRQ( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMRQ( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMRQ( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMRQ( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMRQ( 'L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMRQ( 'R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORMRQ( 'L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMRQ( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMRQ( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMRQ', INFOT, NOUT, LERR, OK )
*
*     DORMR2
*
      SRNAMT = 'DORMR2'
      INFOT = 1
      CALL DORMR2( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORMR2( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORMR2( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORMR2( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMR2( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMR2( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMR2( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMR2( 'L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMR2( 'R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORMR2( 'L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORMR2', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRRQ
*
      END
