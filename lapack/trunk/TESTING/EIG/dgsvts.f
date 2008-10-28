      SUBROUTINE DGSVTS( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
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
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), ALPHA( * ),
     $                   B( LDB, * ), BETA( * ), BF( LDB, * ),
     $                   Q( LDQ, * ), R( LDR, * ), RESULT( 6 ),
     $                   RWORK( * ), U( LDU, * ), V( LDV, * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DGSVTS tests DGGSVD, which computes the GSVD of an M-by-N matrix A
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
*          The M-by-N matrix A.
*
*  AF      (output) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the GSVD of A and B, as returned by DGGSVD,
*          see DGGSVD for further details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*          LDA >= max( 1,M ).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,P)
*          On entry, the P-by-N matrix B.
*
*  BF      (output) DOUBLE PRECISION array, dimension (LDB,N)
*          Details of the GSVD of A and B, as returned by DGGSVD,
*          see DGGSVD for further details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B and BF.
*          LDB >= max(1,P).
*
*  U       (output) DOUBLE PRECISION array, dimension(LDU,M)
*          The M by M orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M).
*
*  V       (output) DOUBLE PRECISION array, dimension(LDV,M)
*          The P by P orthogonal matrix V.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P).
*
*  Q       (output) DOUBLE PRECISION array, dimension(LDQ,N)
*          The N by N orthogonal matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          The generalized singular value pairs of A and B, the
*          ``diagonal'' matrices D1 and D2 are constructed from
*          ALPHA and BETA, see subroutine DGGSVD for details.
*
*  R       (output) DOUBLE PRECISION array, dimension(LDQ,N)
*          The upper triangular matrix R.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,N).
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK,
*          LWORK >= max(M,P,N)*max(M,P,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(M,P,N))
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (6)
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
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K, L
      DOUBLE PRECISION   ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGGSVD, DLACPY, DLASET, DSYRK
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
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL DLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( DLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL DGGSVD( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB,
     $             ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK,
     $             INFO )
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
      CALL DGEMM( 'No transpose', 'No transpose', M, N, N, ONE, A, LDA,
     $            Q, LDQ, ZERO, WORK, LDA )
*
      CALL DGEMM( 'Transpose', 'No transpose', M, N, M, ONE, U, LDU,
     $            WORK, LDA, ZERO, A, LDA )
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
      RESID = DLANGE( '1', M, N, A, LDA, RWORK )
*
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) /
     $                 ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute B := V'*B*Q - D2*R
*
      CALL DGEMM( 'No transpose', 'No transpose', P, N, N, ONE, B, LDB,
     $            Q, LDQ, ZERO, WORK, LDB )
*
      CALL DGEMM( 'Transpose', 'No transpose', P, N, P, ONE, V, LDV,
     $            WORK, LDB, ZERO, B, LDB )
*
      DO 100 I = 1, L
         DO 90 J = I, L
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
   90    CONTINUE
  100 CONTINUE
*
*     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
*
      RESID = DLANGE( '1', P, N, B, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) /
     $                 ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - U'*U
*
      CALL DLASET( 'Full', M, M, ZERO, ONE, WORK, LDQ )
      CALL DSYRK( 'Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK,
     $            LDU )
*
*     Compute norm( I - U'*U ) / ( M * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP
*
*     Compute I - V'*V
*
      CALL DLASET( 'Full', P, P, ZERO, ONE, WORK, LDV )
      CALL DSYRK( 'Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK,
     $            LDV )
*
*     Compute norm( I - V'*V ) / ( P * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
*     Compute I - Q'*Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, WORK, LDQ )
      CALL DSYRK( 'Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK,
     $            LDQ )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
*
*     Check sorting
*
      CALL DCOPY( N, ALPHA, 1, WORK, 1 )
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         IF( I.NE.J ) THEN
            TEMP = WORK( I )
            WORK( I ) = WORK( J )
            WORK( J ) = TEMP
         END IF
  110 CONTINUE
*
      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( WORK( I ).LT.WORK( I+1 ) )
     $      RESULT( 6 ) = ULPINV
  120 CONTINUE
*
      RETURN
*
*     End of DGSVTS
*
      END
