*> \brief \b SERRKT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SERRKT( PATH, NUNIT )
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
*> SERRKT tests the error exits for SKYTRD, SKTEQR and SKYEV.
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SERRKT( PATH, NUNIT )
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
*     NMAX has to be at least 3 or LIW may be too small
*     .. Parameters ..
      INTEGER            NMAX, LIW, LW
      PARAMETER          ( NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J, M, N, NSPLIT, NT
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), Q( NMAX, NMAX ), R( NMAX ),
     $                   TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, SKTEQR, SKYEV, SKTEV, SKYTRD,
     $                   SKTEBZ, SKTEIN, SKYEVX, SKTEVX
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
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, NMAX
         D( J ) = REAL( J )
         E( J ) = 0.0
         I1( J ) = J
         I2( J ) = J
         TAU( J ) = 1.
   30 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the KT path.
*
      IF( LSAMEN( 2, C2, 'KT' ) ) THEN
*
*        SKYTRD
*
         SRNAMT = 'SKYTRD'
         INFOT = 1
         CALL SKYTRD( '/', 0, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRD( 'U', -1, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTRD( 'U', 2, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKYTRD( 'U', 0, A, 1, E, TAU, W, 0, INFO )
         CALL CHKXER( 'SKYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SKTEBZ
*
         SRNAMT = 'SKTEBZ'
         INFOT = 1
         CALL SKTEBZ( '/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKTEBZ( 'A', '/', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKTEBZ( 'A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SKTEBZ( 'V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKTEBZ( 'I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKTEBZ( 'I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SKTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        SKTEIN
*
         SRNAMT = 'SKTEIN'
         INFOT = 1
         CALL SKTEIN( -1, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKTEIN( 0, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKTEIN( 0, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKTEIN( 2, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SKTEQR
*
         SRNAMT = 'SKTEQR'
         INFOT = 1
         CALL SKTEQR( '/', 0, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKTEQR( 'N', -1, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SKTEQR( 'V', 2, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SKYEV
*
         SRNAMT = 'SKYEV '
         INFOT = 1
         CALL SKYEV( '/', 'U', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYEV( 'N', '/', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKYEV( 'N', 'U', -1, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SKYEV( 'N', 'U', 2, A, 1, X, W, 3, INFO )
         CALL CHKXER( 'SKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKYEV( 'N', 'U', 2, A, 2, X, W, 2, INFO )
         CALL CHKXER( 'SKYEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 5
*
*        SKTEV
*
         SRNAMT = 'SKTEV '
         INFOT = 1
         CALL SKTEV( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKTEV( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKTEV( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SKTEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SKYEVX
*
         SRNAMT = 'SKYEVX'
         INFOT = 1
         CALL SKYEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYEVX( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKYEVX( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         INFOT = 4
         CALL SKYEVX( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M,
     $                X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKYEVX( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKYEVX( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKYEVX( 'N', 'V', 'U', 1, A, 1, -2.0, -1.0, 0, 0, 0.0, M,
     $                X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SKYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SKYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SKYEVX( 'N', 'I', 'U', 3, A, 3, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SKYEVX( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 1, 2, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL SKYEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL SKYEVX( 'V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 0, IW, I3, INFO )
         CALL CHKXER( 'SKYEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        SKTEVX
*
         SRNAMT = 'SKTEVX'
         INFOT = 1
         CALL SKTEVX( '/', 'A', 0, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKTEVX( 'N', '/', 0, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKTEVX( 'N', 'A', -1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKTEVX( 'N', 'V', 1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SKTEVX( 'N', 'V', 1, E, -2.0, -1.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKTEVX( 'N', 'I', 1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKTEVX( 'N', 'I', 1, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKTEVX( 'N', 'I', 3, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKTEVX( 'N', 'I', 2, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL SKTEVX( 'V', 'A', 2, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'SKTEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 10
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits',
     $      ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
*
      RETURN
*
*     End of SERRKT
*
      END
