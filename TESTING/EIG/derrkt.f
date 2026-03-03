*> \brief \b DERRKT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRKT( PATH, NUNIT )
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
*> DERRKT tests the error exits for DKYTRD, DKTEQR and DKYEV.
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DERRKT( PATH, NUNIT )
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
      DOUBLE PRECISION   A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), Q( NMAX, NMAX ), R( NMAX ),
     $                   TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, DKTEQR, DKYEV, DKTEV, DKYTRD,
     $                   DKTEBZ, DKTEIN, DKYEVX, DKTEVX
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
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / DBLE( I+J )
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, NMAX
         D( J ) = DBLE( J )
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
*        DKYTRD
*
         SRNAMT = 'DKYTRD'
         INFOT = 1
         CALL DKYTRD( '/', 0, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRD( 'U', -1, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTRD( 'U', 2, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DKYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKYTRD( 'U', 0, A, 1, E, TAU, W, 0, INFO )
         CALL CHKXER( 'DKYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DKTEBZ
*
         SRNAMT = 'DKTEBZ'
         INFOT = 1
         CALL DKTEBZ( '/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKTEBZ( 'A', '/', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKTEBZ( 'A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DKTEBZ( 'V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKTEBZ( 'I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKTEBZ( 'I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DKTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        DKTEIN
*
         SRNAMT = 'DKTEIN'
         INFOT = 1
         CALL DKTEIN( -1, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKTEIN( 0, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKTEIN( 0, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKTEIN( 2, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DKTEQR
*
         SRNAMT = 'DKTEQR'
         INFOT = 1
         CALL DKTEQR( '/', 0, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKTEQR( 'N', -1, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DKTEQR( 'V', 2, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DKYEV
*
         SRNAMT = 'DKYEV '
         INFOT = 1
         CALL DKYEV( '/', 'U', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYEV( 'N', '/', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKYEV( 'N', 'U', -1, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DKYEV( 'N', 'U', 2, A, 1, X, W, 3, INFO )
         CALL CHKXER( 'DKYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKYEV( 'N', 'U', 2, A, 2, X, W, 2, INFO )
         CALL CHKXER( 'DKYEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 5
*
*        DKTEV
*
         SRNAMT = 'DKTEV '
         INFOT = 1
         CALL DKTEV( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKTEV( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKTEV( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DKTEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DKYEVX
*
         SRNAMT = 'DKYEVX'
         INFOT = 1
         CALL DKYEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYEVX( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKYEVX( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 1, IW, I3, INFO )
         INFOT = 4
         CALL DKYEVX( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M,
     $                X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKYEVX( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKYEVX( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKYEVX( 'N', 'V', 'U', 1, A, 1, -2.0, -1.0, 0, 0, 0.0, M,
     $                X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DKYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DKYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DKYEVX( 'N', 'I', 'U', 3, A, 3, 0.0, 0.0, 2, 1, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DKYEVX( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 1, 2, 0.0, M, X,
     $                Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DKYEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DKYEVX( 'V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
     $                Z, 1, W, 0, IW, I3, INFO )
         CALL CHKXER( 'DKYEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        DKTEVX
*
         SRNAMT = 'DKTEVX'
         INFOT = 1
         CALL DKTEVX( '/', 'A', 0, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKTEVX( 'N', '/', 0, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKTEVX( 'N', 'A', -1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKTEVX( 'N', 'V', 1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DKTEVX( 'N', 'V', 1, E, -2.0, -1.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKTEVX( 'N', 'I', 1, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKTEVX( 'N', 'I', 1, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKTEVX( 'N', 'I', 3, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKTEVX( 'N', 'I', 2, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL DKTEVX( 'V', 'A', 2, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z,
     $                1, W, IW, I3, INFO )
         CALL CHKXER( 'DKTEVX', INFOT, NOUT, LERR, OK )
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
*     End of DERRKT
*
      END
