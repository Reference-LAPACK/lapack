      SUBROUTINE CCSDTS( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T,
     $                   LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK,
     $                   RWORK, RESULT )
      IMPLICIT NONE
*
*     Originally xGSVTS
*  -- LAPACK test routine (version 3.3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     Adapted to CCSDTS by
*     July 2010
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               RESULT( 9 ), RWORK( * ), THETA( * )
      COMPLEX            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
     $                   V2T( LDV2T, * ), WORK( LWORK ), X( LDX, * ),
     $                   XF( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  CCSDTS tests CUNCSD, which, given an M-by-M partitioned unitary
*  matrix X,
*               Q  M-Q
*        X = [ X11 X12 ] P   ,
*            [ X21 X22 ] M-P
*
*  computes the CSD
*
*        [ U1    ]**T * [ X11 X12 ] * [ V1    ]
*        [    U2 ]      [ X21 X22 ]   [    V2 ]
*
*                              [  I  0  0 |  0  0  0 ]
*                              [  0  C  0 |  0 -S  0 ]
*                              [  0  0  0 |  0  0 -I ]
*                            = [---------------------] = [ D11 D12 ] .
*                              [  0  0  0 |  I  0  0 ]   [ D21 D22 ]
*                              [  0  S  0 |  0  C  0 ]
*                              [  0  0  I |  0  0  0 ]
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix X.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix X11.  P >= 0.
*
*  Q       (input) INTEGER
*          The number of columns of the matrix X11.  Q >= 0.
*
*  X       (input) COMPLEX array, dimension (LDX,M)
*          The M-by-M matrix X.
*
*  XF      (output) COMPLEX array, dimension (LDX,M)
*          Details of the CSD of X, as returned by CUNCSD;
*          see CUNCSD for further details.
*
*  LDX     (input) INTEGER
*          The leading dimension of the arrays X and XF.
*          LDX >= max( 1,M ).
*
*  U1      (output) COMPLEX array, dimension(LDU1,P)
*          The P-by-P unitary matrix U1.
*
*  LDU1    (input) INTEGER
*          The leading dimension of the array U1. LDU >= max(1,P).
*
*  U2      (output) COMPLEX array, dimension(LDU2,M-P)
*          The (M-P)-by-(M-P) unitary matrix U2.
*
*  LDU2    (input) INTEGER
*          The leading dimension of the array U2. LDU >= max(1,M-P).
*
*  V1T     (output) COMPLEX array, dimension(LDV1T,Q)
*          The Q-by-Q unitary matrix V1T.
*
*  LDV1T   (input) INTEGER
*          The leading dimension of the array V1T. LDV1T >=
*          max(1,Q).
*
*  V2T     (output) COMPLEX array, dimension(LDV2T,M-Q)
*          The (M-Q)-by-(M-Q) unitary matrix V2T.
*
*  LDV2T   (input) INTEGER
*          The leading dimension of the array V2T. LDV2T >=
*          max(1,M-Q).
*
*  THETA   (output) REAL array, dimension MIN(P,M-P,Q,M-Q)
*          The CS values of X; the essentially diagonal matrices C and
*          S are constructed from THETA; see subroutine CUNCSD for
*          details.
*
*  IWORK   (workspace) INTEGER array, dimension (M)
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK
*
*  RWORK   (workspace) REAL array
*
*  RESULT  (output) REAL array, dimension (9)
*          The test ratios:
*          RESULT(1) = norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 )
*          RESULT(2) = norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 )
*          RESULT(3) = norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 )
*          RESULT(4) = norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 )
*          RESULT(5) = norm( I - U1'*U1 ) / ( MAX(1,P)*ULP )
*          RESULT(6) = norm( I - U2'*U2 ) / ( MAX(1,M-P)*ULP )
*          RESULT(7) = norm( I - V1T'*V1T ) / ( MAX(1,Q)*ULP )
*          RESULT(8) = norm( I - V2T'*V2T ) / ( MAX(1,M-Q)*ULP )
*          RESULT(9) = 0        if THETA is in increasing order and
*                               all angles are in [0,pi/2];
*                    = ULPINV   otherwise.
*          ( EPS2 = MAX( norm( I - X'*X ) / M, ULP ). )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               PIOVER2, REALONE, REALZERO
      PARAMETER          ( PIOVER2 = 1.57079632679489662E0,
     $                     REALONE = 1.0E0, REALZERO = 0.0E0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = (0.0E0,0.0E0), ONE = (1.0E0,0.0E0) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, R
      REAL               EPS2, RESID, ULP, ULPINV
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CLANGE, CLANHE
      EXTERNAL           SLAMCH, CLANGE, CLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLASET, CUNCSD, CHERK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      ULPINV = REALONE / ULP
      CALL CLASET( 'Full', M, M, ZERO, ONE, WORK, LDX )
      CALL CHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, X, LDX,
     $            ONE, WORK, LDX )
      EPS2 = MAX( ULP, 
     $            CLANGE( '1', M, M, WORK, LDX, RWORK ) / REAL( M ) )
      R = MIN( P, M-P, Q, M-Q )
*
*     Copy the matrix X to the array XF.
*
      CALL CLACPY( 'Full', M, M, X, LDX, XF, LDX )
*
*     Compute the CSD
*
      CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'D', M, P, Q, XF(1,1), LDX,
     $             XF(1,Q+1), LDX, XF(P+1,1), LDX, XF(P+1,Q+1), LDX,
     $             THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T,
     $             WORK, LWORK, RWORK, 17*(R+2), IWORK, INFO )
*
*     Compute X := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', P, Q, Q, ONE,
     $            X, LDX, V1T, LDV1T, ZERO, WORK, LDX )
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', P, Q, P, ONE,
     $            U1, LDU1, WORK, LDX, ZERO, X, LDX )
*
      DO I = 1, MIN(P,Q)-R
         X(I,I) = X(I,I) - ONE
      END DO
      DO I = 1, R
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) =
     $           X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) - CMPLX( COS(THETA(I)),
     $              0.0E0 )
      END DO
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', P, M-Q, M-Q,
     $            ONE, X(1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX )
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', P, M-Q, P,
     $            ONE, U1, LDU1, WORK, LDX, ZERO, X(1,Q+1), LDX )
*
      DO I = 1, MIN(P,M-Q)-R
         X(P-I+1,M-I+1) = X(P-I+1,M-I+1) + ONE
      END DO
      DO I = 1, R
         X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) =
     $      X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) +
     $      CMPLX( SIN(THETA(R-I+1)), 0.0E0 )
      END DO
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', M-P, Q, Q, ONE,
     $            X(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX )
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', M-P, Q, M-P,
     $            ONE, U2, LDU2, WORK, LDX, ZERO, X(P+1,1), LDX )
*
      DO I = 1, MIN(M-P,Q)-R
         X(M-I+1,Q-I+1) = X(M-I+1,Q-I+1) - ONE
      END DO
      DO I = 1, R
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) =
     $             X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) -
     $             CMPLX( SIN(THETA(R-I+1)), 0.0E0 )
      END DO
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', M-P, M-Q, M-Q,
     $            ONE, X(P+1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX )
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', M-P, M-Q, M-P,
     $            ONE, U2, LDU2, WORK, LDX, ZERO, X(P+1,Q+1), LDX )
*
      DO I = 1, MIN(M-P,M-Q)-R
         X(P+I,Q+I) = X(P+I,Q+I) - ONE
      END DO
      DO I = 1, R
         X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) =
     $      X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) -
     $      CMPLX( COS(THETA(I)), 0.0E0 )
      END DO
*
*     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
*
      RESID = CLANGE( '1', P, Q, X, LDX, RWORK )
      RESULT( 1 ) = ( RESID / REAL(MAX(1,P,Q)) ) / EPS2
*
*     Compute norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 ) .
*
      RESID = CLANGE( '1', P, M-Q, X(1,Q+1), LDX, RWORK )
      RESULT( 2 ) = ( RESID / REAL(MAX(1,P,M-Q)) ) / EPS2
*
*     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
*
      RESID = CLANGE( '1', M-P, Q, X(P+1,1), LDX, RWORK )
      RESULT( 3 ) = ( RESID / REAL(MAX(1,M-P,Q)) ) / EPS2
*
*     Compute norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 ) .
*
      RESID = CLANGE( '1', M-P, M-Q, X(P+1,Q+1), LDX, RWORK )
      RESULT( 4 ) = ( RESID / REAL(MAX(1,M-P,M-Q)) ) / EPS2
*
*     Compute I - U1'*U1
*
      CALL CLASET( 'Full', P, P, ZERO, ONE, WORK, LDU1 )
      CALL CHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, U1, LDU1,
     $            ONE, WORK, LDU1 )
*
*     Compute norm( I - U'*U ) / ( MAX(1,P) * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', P, WORK, LDU1, RWORK )
      RESULT( 5 ) = ( RESID / REAL(MAX(1,P)) ) / ULP
*
*     Compute I - U2'*U2
*
      CALL CLASET( 'Full', M-P, M-P, ZERO, ONE, WORK, LDU2 )
      CALL CHERK( 'Upper', 'Conjugate transpose', M-P, M-P, -ONE, U2,
     $            LDU2, ONE, WORK, LDU2 )
*
*     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', M-P, WORK, LDU2, RWORK )
      RESULT( 6 ) = ( RESID / REAL(MAX(1,M-P)) ) / ULP
*
*     Compute I - V1T*V1T'
*
      CALL CLASET( 'Full', Q, Q, ZERO, ONE, WORK, LDV1T )
      CALL CHERK( 'Upper', 'No transpose', Q, Q, -ONE, V1T, LDV1T, ONE,
     $            WORK, LDV1T )
*
*     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', Q, WORK, LDV1T, RWORK )
      RESULT( 7 ) = ( RESID / REAL(MAX(1,Q)) ) / ULP
*
*     Compute I - V2T*V2T'
*
      CALL CLASET( 'Full', M-Q, M-Q, ZERO, ONE, WORK, LDV2T )
      CALL CHERK( 'Upper', 'No transpose', M-Q, M-Q, -ONE, V2T, LDV2T,
     $            ONE, WORK, LDV2T )
*
*     Compute norm( I - V2T*V2T' ) / ( MAX(1,M-Q) * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', M-Q, WORK, LDV2T, RWORK )
      RESULT( 8 ) = ( RESID / REAL(MAX(1,M-Q)) ) / ULP
*
*     Check sorting
*
      RESULT(9) = REALZERO
      DO I = 1, R
         IF( THETA(I).LT.REALZERO .OR. THETA(I).GT.PIOVER2 ) THEN
            RESULT(9) = ULPINV
         END IF
         IF( I.GT.1) THEN
            IF ( THETA(I).LT.THETA(I-1) ) THEN
               RESULT(9) = ULPINV
            END IF
         END IF
      END DO
*
      RETURN
*      
*     End of CCSDTS
*
      END

