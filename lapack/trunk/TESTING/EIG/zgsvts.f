      SUBROUTINE ZGSVTS( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
     $                   LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
     $                   LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), Q( LDQ, * ), R( LDR, * ),
     $                   U( LDU, * ), V( LDV, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZGSVTS tests ZGGSVD, which computes the GSVD of an M-by-N matrix A
*  and a P-by-N matrix B:
*               U'*A*Q = D1*R and V'*B*Q = D2*R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,M)
*          The M-by-N matrix A.
*
*  AF      (output) COMPLEX*16 array, dimension (LDA,N)
*          Details of the GSVD of A and B, as returned by ZGGSVD,
*          see ZGGSVD for further details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*          LDA >= max( 1,M ).
*
*  B       (input) COMPLEX*16 array, dimension (LDB,P)
*          On entry, the P-by-N matrix B.
*
*  BF      (output) COMPLEX*16 array, dimension (LDB,N)
*          Details of the GSVD of A and B, as returned by ZGGSVD,
*          see ZGGSVD for further details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B and BF.
*          LDB >= max(1,P).
*
*  U       (output) COMPLEX*16 array, dimension(LDU,M)
*          The M by M unitary matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M).
*
*  V       (output) COMPLEX*16 array, dimension(LDV,M)
*          The P by P unitary matrix V.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P).
*
*  Q       (output) COMPLEX*16 array, dimension(LDQ,N)
*          The N by N unitary matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          The generalized singular value pairs of A and B, the
*          ``diagonal'' matrices D1 and D2 are constructed from
*          ALPHA and BETA, see subroutine ZGGSVD for details.
*
*  R       (output) COMPLEX*16 array, dimension(LDQ,N)
*          The upper triangular matrix R.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,N).
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK,
*          LWORK >= max(M,P,N)*max(M,P,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(M,P,N))
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (5)
*          The test ratios:
*          RESULT(1) = norm( U'*A*Q - D1*R ) / ( MAX(M,N)*norm(A)*ULP)
*          RESULT(2) = norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP)
*          RESULT(3) = norm( I - U'*U ) / ( M*ULP )
*          RESULT(4) = norm( I - V'*V ) / ( P*ULP )
*          RESULT(5) = norm( I - Q'*Q ) / ( N*ULP )
*          RESULT(6) = 0        if ALPHA is in decreasing order;
*                    = ULPINV   otherwise.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K, L
      DOUBLE PRECISION   ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
      EXTERNAL           DLAMCH, ZLANGE, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, ZGEMM, ZGGSVD, ZHERK, ZLACPY, ZLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = DLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL ZLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL ZGGSVD( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB,
     $             ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, RWORK,
     $             IWORK, INFO )
*
*     Copy R
*
      DO 20 I = 1, MIN( K+L, M )
         DO 10 J = I, K + L
            R( I, J ) = AF( I, N-K-L+J )
   10    CONTINUE
   20 CONTINUE
*
      IF( M-K-L.LT.0 ) THEN
         DO 40 I = M + 1, K + L
            DO 30 J = I, K + L
               R( I, J ) = BF( I-K, N-K-L+J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Compute A:= U'*A*Q - D1*R
*
      CALL ZGEMM( 'No transpose', 'No transpose', M, N, N, CONE, A, LDA,
     $            Q, LDQ, CZERO, WORK, LDA )
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M, CONE,
     $            U, LDU, WORK, LDA, CZERO, A, LDA )
*
      DO 60 I = 1, K
         DO 50 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
   50    CONTINUE
   60 CONTINUE
*
      DO 80 I = K + 1, MIN( K+L, M )
         DO 70 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - ALPHA( I )*R( I, J )
   70    CONTINUE
   80 CONTINUE
*
*     Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
*
      RESID = ZLANGE( '1', M, N, A, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) /
     $                 ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute B := V'*B*Q - D2*R
*
      CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, CONE, B, LDB,
     $            Q, LDQ, CZERO, WORK, LDB )
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', P, N, P, CONE,
     $            V, LDV, WORK, LDB, CZERO, B, LDB )
*
      DO 100 I = 1, L
         DO 90 J = I, L
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
   90    CONTINUE
  100 CONTINUE
*
*     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
*
      RESID = ZLANGE( '1', P, N, B, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) /
     $                 ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - U'*U
*
      CALL ZLASET( 'Full', M, M, CZERO, CONE, WORK, LDQ )
      CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, U, LDU,
     $            ONE, WORK, LDU )
*
*     Compute norm( I - U'*U ) / ( M * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP
*
*     Compute I - V'*V
*
      CALL ZLASET( 'Full', P, P, CZERO, CONE, WORK, LDV )
      CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, V, LDV,
     $            ONE, WORK, LDV )
*
*     Compute norm( I - V'*V ) / ( P * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
*     Compute I - Q'*Q
*
      CALL ZLASET( 'Full', N, N, CZERO, CONE, WORK, LDQ )
      CALL ZHERK( 'Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ,
     $            ONE, WORK, LDQ )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
*
*     Check sorting
*
      CALL DCOPY( N, ALPHA, 1, RWORK, 1 )
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         IF( I.NE.J ) THEN
            TEMP = RWORK( I )
            RWORK( I ) = RWORK( J )
            RWORK( J ) = TEMP
         END IF
  110 CONTINUE
*
      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( RWORK( I ).LT.RWORK( I+1 ) )
     $      RESULT( 6 ) = ULPINV
  120 CONTINUE
*
      RETURN
*
*     End of ZGSVTS
*
      END
