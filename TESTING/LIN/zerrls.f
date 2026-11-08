*> \brief \b ZERRLS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZERRLS( PATH, NUNIT )
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
*> ZERRLS tests the error exits for the COMPLEX*16 least squares
*> driver routines (ZGELS, ZGELST, ZGETSLS, CGELSS, CGELSY, CGELSD).
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZERRLS( PATH, NUNIT )
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
      INTEGER            NMAX, LW, LIW
      PARAMETER          ( NMAX = 2, LW = 1000, LIW = 30 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      LOGICAL            NANB
      INTEGER            INFO, IRNK, J
      DOUBLE PRECISION   INF, NAN, ONE, RCOND, RZERO, TOL
*     ..
*     .. Local Arrays ..
      INTEGER            IP( LIW )
      DOUBLE PRECISION   E( NMAX )
      DOUBLE PRECISION   RW( LW ), S( NMAX )
      COMPLEX*16         A( NMAX, NMAX ), B( NMAX, NMAX ), W( LW )
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAMEN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZGELS, ZGELSD, ZGELSS, ZGELST,
     $                   ZGELSY, ZGETSLS, ZLALSD
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
      INTRINSIC          ABS, DBLE, DIMAG, LEN_TRIM, SQRT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
      A( 1, 2 ) = ( 2.0D+0, 0.0D+0 )
      A( 2, 2 ) = ( 3.0D+0, 0.0D+0 )
      A( 2, 1 ) = ( 4.0D+0, 0.0D+0 )
      OK = .TRUE.
      WRITE( NOUT, FMT = * )
*
*     Test error exits for the least squares driver routines.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        ZGELS
*
         SRNAMT = 'ZGELS '
         INFOT = 1
         CALL ZGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELS', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
*
*        ZGELST
*
         SRNAMT = 'ZGELST'
         INFOT = 1
         CALL ZGELST( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELST( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELST( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGELST( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGELST( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELST( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELST( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGELST( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
*
*        ZGETSLS
*
         SRNAMT = 'ZGETSLS'
         INFOT = 1
         CALL ZGETSLS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGETSLS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGETSLS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGETSLS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGETSLS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGETSLS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGETSLS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
*
*        ZGELSS
*
         SRNAMT = 'ZGELSS'
         INFOT = 1
         CALL ZGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
*
*        ZGELSY
*
         SRNAMT = 'ZGELSY'
         INFOT = 1
         CALL ZGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZGELSY( 0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
*
*        ZGELSD
*
         SRNAMT = 'ZGELSD'
         INFOT = 1
         CALL ZGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW,
     $                IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW,
     $                IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW,
     $                IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW,
     $                IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW,
     $                IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP,
     $                INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
*
*        ZGELSD must return INFO = -4, without calling XERBLA, when A
*        contains a NaN or an infinity.
*
         ONE = 1.0D+0
         RZERO = 0.0D+0
         NAN = SQRT( -ONE )
         INF = ONE / RZERO
         RCOND = -1.0D+0
         DO 10 J = 1, 2
            A( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
            A( 1, 2 ) = ( 1.0D+0, 0.0D+0 )
            A( 2, 1 ) = ( 1.0D+0, 0.0D+0 )
            A( 2, 2 ) = NAN
            IF( J.EQ.2 )
     $         A( 2, 2 ) = INF
            B( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
            B( 2, 1 ) = ( 1.0D+0, 0.0D+0 )
            SRNAMT = 'ZGELSD'
            INFOT = 0
            LERR = .FALSE.
            CALL ZGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, LW, RW,
     $                   IP, INFO )
            IF( LERR .OR. INFO.NE.-4 ) THEN
               WRITE( NOUT, FMT = 9999 )'ZGELSD', -4, INFO
               OK = .FALSE.
            END IF
   10    CONTINUE
*
*        ZLALSD must return INFO = -5 or -6, without calling XERBLA,
*        when D or E contains a NaN or an infinity.
*
         DO 20 J = 1, 2
            S( 1 ) = ONE
            S( 2 ) = ONE
            E( 1 ) = ONE
            IF( J.EQ.1 ) THEN
               S( 2 ) = INF
            ELSE
               E( 1 ) = NAN
            END IF
            B( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
            B( 2, 1 ) = ( 1.0D+0, 0.0D+0 )
            SRNAMT = 'ZLALSD'
            INFOT = 0
            LERR = .FALSE.
            CALL ZLALSD( 'U', 25, 2, 1, S, E, B, 2, RCOND, IRNK, W, RW,
     $                   IP, INFO )
            IF( LERR .OR. INFO.NE.-4-J ) THEN
               WRITE( NOUT, FMT = 9999 )'ZLALSD', -4-J, INFO
               OK = .FALSE.
            END IF
   20    CONTINUE
*
*        A NaN in one column of B must leave the solution for the
*        other column intact: each driver returns INFO = 0, without
*        calling XERBLA, the solution for the finite column, and a
*        NaN in the solution for the other one.  The drivers do not
*        check B for a NaN or an infinity.
*
         TOL = 100*DLAMCH( 'Epsilon' )
         DO 30 J = 1, 6
            A( 1, 1 ) = ( 3.0D+0, 0.0D+0 )
            A( 2, 1 ) = ( 1.0D+0, 0.0D+0 )
            A( 1, 2 ) = ( 1.0D+0, 0.0D+0 )
            A( 2, 2 ) = ( 2.0D+0, 0.0D+0 )
            B( 1, 1 ) = ( 4.0D+0, 0.0D+0 )
            B( 2, 1 ) = ( 3.0D+0, 0.0D+0 )
            B( 1, 2 ) = ( 1.0D+0, 0.0D+0 )
            B( 2, 2 ) = NAN
            IP( 1 ) = 0
            IP( 2 ) = 0
            INFOT = 0
            LERR = .FALSE.
            IF( J.EQ.1 ) THEN
               SRNAMT = 'ZGELS'
               CALL ZGELS( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.2 ) THEN
               SRNAMT = 'ZGELST'
               CALL ZGELST( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.3 ) THEN
               SRNAMT = 'ZGETSLS'
               CALL ZGETSLS( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.4 ) THEN
               SRNAMT = 'ZGELSS'
               CALL ZGELSS( 2, 2, 2, A, 2, B, 2, S, RCOND, IRNK, W, LW,
     $                      RW, INFO )
            ELSE IF( J.EQ.5 ) THEN
               SRNAMT = 'ZGELSY'
               CALL ZGELSY( 2, 2, 2, A, 2, B, 2, IP, RCOND, IRNK, W, LW,
     $                      RW, INFO )
            ELSE
               SRNAMT = 'ZGELSD'
               CALL ZGELSD( 2, 2, 2, A, 2, B, 2, S, RCOND, IRNK, W, LW,
     $                      RW, IP, INFO )
            END IF
            NANB = DISNAN( DBLE( B( 1, 2 ) ) ) .OR.
     $             DISNAN( DIMAG( B( 1, 2 ) ) ) .OR.
     $             DISNAN( DBLE( B( 2, 2 ) ) ) .OR.
     $             DISNAN( DIMAG( B( 2, 2 ) ) )
            IF( LERR .OR. INFO.NE.0 .OR. .NOT.NANB .OR.
     $          .NOT.( ABS( B( 1, 1 )-ONE ).LE.TOL .AND.
     $                 ABS( B( 2, 1 )-ONE ).LE.TOL ) ) THEN
               WRITE( NOUT, FMT = 9998 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ),
     $            INFO
               OK = .FALSE.
            END IF
   30    CONTINUE
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
 9999 FORMAT( ' *** ', A6, ' did not return INFO = ', I2, ' without',
     $      ' calling XERBLA for input with a NaN or an infinity',
     $      ' (INFO = ', I6, ') ***' )
 9998 FORMAT( ' *** ', A, ' with a NaN in one column of B did not',
     $      ' return INFO = 0 and the solution for the other column',
     $      ' without calling XERBLA (INFO = ', I6, ') ***' )
*
*     End of ZERRLS
*
      END
