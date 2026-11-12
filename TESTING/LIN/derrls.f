*> \brief \b DERRLS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRLS( PATH, NUNIT )
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
*> DERRLS tests the error exits for the DOUBLE PRECISION least squares
*> driver routines (DGELS, DGELST, DGETSLS, SGELSS, SGELSY, SGELSD).
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DERRLS( PATH, NUNIT )
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
      DOUBLE PRECISION   A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ),
     $                   W( LW )
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAMEN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELST,
     $                   DGELSY, DGETSLS, DLALSD
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
      INTRINSIC          ABS, LEN_TRIM, SQRT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.0D+0
      A( 1, 2 ) = 2.0D+0
      A( 2, 2 ) = 3.0D+0
      A( 2, 1 ) = 4.0D+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        Test error exits for the least squares driver routines.
*
*        DGELS
*
         SRNAMT = 'DGELS '
         INFOT = 1
         CALL DGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGELS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGELS', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
*
*        DGELST
*
         SRNAMT = 'DGELST'
         INFOT = 1
         CALL DGELST( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELST( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELST( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGELST( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGELST( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGELST( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGELST( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGELST( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELST', INFOT, NOUT, LERR, OK )
*
*        DGETSLS
*
         SRNAMT = 'DGETSLS'
         INFOT = 1
         CALL DGETSLS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGETSLS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGETSLS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGETSLS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGETSLS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGETSLS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGETSLS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGETSLS', INFOT, NOUT, LERR, OK )
*
*        DGELSS
*
         SRNAMT = 'DGELSS'
         INFOT = 1
         CALL DGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
*
*        DGELSY
*
         SRNAMT = 'DGELSY'
         INFOT = 1
         CALL DGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGELSY( 2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
*
*        DGELSD
*
         SRNAMT = 'DGELSD'
         INFOT = 1
         CALL DGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
*
*        DGELSD must return INFO = -4, without calling XERBLA, when A
*        contains a NaN or an infinity.
*
         ONE = 1.0D+0
         RZERO = 0.0D+0
         NAN = SQRT( -ONE )
         INF = ONE / RZERO
         RCOND = -1.0D+0
         DO 10 J = 1, 2
            A( 1, 1 ) = 1.0D+0
            A( 1, 2 ) = 1.0D+0
            A( 2, 1 ) = 1.0D+0
            A( 2, 2 ) = NAN
            IF( J.EQ.2 )
     $         A( 2, 2 ) = INF
            B( 1, 1 ) = 1.0D+0
            B( 2, 1 ) = 1.0D+0
            SRNAMT = 'DGELSD'
            INFOT = 0
            LERR = .FALSE.
            CALL DGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, LW,
     $                   IP, INFO )
            IF( LERR .OR. INFO.NE.-4 ) THEN
               WRITE( NOUT, FMT = 9999 )'DGELSD', -4, INFO
               OK = .FALSE.
            END IF
   10    CONTINUE
*
*        DLALSD must return INFO = -5 or -6, without calling XERBLA,
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
            B( 1, 1 ) = 1.0D+0
            B( 2, 1 ) = 1.0D+0
            SRNAMT = 'DLALSD'
            INFOT = 0
            LERR = .FALSE.
            CALL DLALSD( 'U', 25, 2, 1, S, E, B, 2, RCOND, IRNK, W,
     $                   IP, INFO )
            IF( LERR .OR. INFO.NE.-4-J ) THEN
               WRITE( NOUT, FMT = 9999 )'DLALSD', -4-J, INFO
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
            A( 1, 1 ) = 3.0D+0
            A( 2, 1 ) = 1.0D+0
            A( 1, 2 ) = 1.0D+0
            A( 2, 2 ) = 2.0D+0
            B( 1, 1 ) = 4.0D+0
            B( 2, 1 ) = 3.0D+0
            B( 1, 2 ) = 1.0D+0
            B( 2, 2 ) = NAN
            IP( 1 ) = 0
            IP( 2 ) = 0
            INFOT = 0
            LERR = .FALSE.
            IF( J.EQ.1 ) THEN
               SRNAMT = 'DGELS'
               CALL DGELS( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.2 ) THEN
               SRNAMT = 'DGELST'
               CALL DGELST( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.3 ) THEN
               SRNAMT = 'DGETSLS'
               CALL DGETSLS( 'N', 2, 2, 2, A, 2, B, 2, W, LW, INFO )
            ELSE IF( J.EQ.4 ) THEN
               SRNAMT = 'DGELSS'
               CALL DGELSS( 2, 2, 2, A, 2, B, 2, S, RCOND, IRNK, W, LW,
     $                      INFO )
            ELSE IF( J.EQ.5 ) THEN
               SRNAMT = 'DGELSY'
               CALL DGELSY( 2, 2, 2, A, 2, B, 2, IP, RCOND, IRNK, W, LW,
     $                      INFO )
            ELSE
               SRNAMT = 'DGELSD'
               CALL DGELSD( 2, 2, 2, A, 2, B, 2, S, RCOND, IRNK, W, LW,
     $                      IP, INFO )
            END IF
            NANB = DISNAN( B( 1, 2 ) ) .OR. DISNAN( B( 2, 2 ) )
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
*     End of DERRLS
*
      END
