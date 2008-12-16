      SUBROUTINE ZGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
     $                   BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RESULT( 4 ), RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ),
     $                   R( LDA, * ), T( LDB, * ), TAUA( * ), TAUB( * ),
     $                   WORK( LWORK ), Z( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGRQTS tests ZGGRQF, which computes the GRQ factorization of an
*  M-by-N matrix A and a P-by-N matrix B: A = R*Q and B = Z*T*Q.
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
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The M-by-N matrix A.
*
*  AF      (output) COMPLEX*16 array, dimension (LDA,N)
*          Details of the GRQ factorization of A and B, as returned
*          by ZGGRQF, see CGGRQF for further details.
*
*  Q       (output) COMPLEX*16 array, dimension (LDA,N)
*          The N-by-N unitary matrix Q.
*
*  R       (workspace) COMPLEX*16 array, dimension (LDA,MAX(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, R and Q.
*          LDA >= max(M,N).
*
*  TAUA    (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by DGGQRC.
*
*  B       (input) COMPLEX*16 array, dimension (LDB,N)
*          On entry, the P-by-N matrix A.
*
*  BF      (output) COMPLEX*16 array, dimension (LDB,N)
*          Details of the GQR factorization of A and B, as returned
*          by ZGGRQF, see CGGRQF for further details.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDB,P)
*          The P-by-P unitary matrix Z.
*
*  T       (workspace) COMPLEX*16 array, dimension (LDB,max(P,N))
*
*  BWK     (workspace) COMPLEX*16 array, dimension (LDB,N)
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B, BF, Z and T.
*          LDB >= max(P,N).
*
*  TAUB    (output) COMPLEX*16 array, dimension (min(P,N))
*          The scalar factors of the elementary reflectors, as returned
*          by DGGRQF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK, LWORK >= max(M,P,N)**2.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (4)
*          The test ratios:
*            RESULT(1) = norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP)
*            RESULT(2) = norm( T*Q - Z'*B ) / (MAX(P,N)*norm(B)*ULP)
*            RESULT(3) = norm( I - Q'*Q ) / ( N*ULP )
*            RESULT(4) = norm( I - Z'*Z ) / ( P*ULP )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         CROGUE
      PARAMETER          ( CROGUE = ( -1.0D+10, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      DOUBLE PRECISION   ANORM, BNORM, RESID, ULP, UNFL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
      EXTERNAL           DLAMCH, ZLANGE, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZGGRQF, ZHERK, ZLACPY, ZLASET, ZUNGQR,
     $                   ZUNGRQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ULP = DLAMCH( 'Precision' )
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
      CALL ZGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK,
     $             INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL ZLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N )
     $      CALL ZLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )
         IF( M.GT.1 )
     $      CALL ZLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA,
     $                   Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 )
     $      CALL ZLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA,
     $                   Q( 2, 1 ), LDA )
      END IF
      CALL ZUNGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL ZLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
      IF( P.GT.1 )
     $   CALL ZLACPY( 'Lower', P-1, N, BF( 2, 1 ), LDB, Z( 2, 1 ), LDB )
      CALL ZUNGQR( P, P, MIN( P, N ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL ZLASET( 'Full', M, N, CZERO, CZERO, R, LDA )
      IF( M.LE.N ) THEN
         CALL ZLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ),
     $                LDA )
      ELSE
         CALL ZLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         CALL ZLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ),
     $                LDA )
      END IF
*
*     Copy T
*
      CALL ZLASET( 'Full', P, N, CZERO, CZERO, T, LDB )
      CALL ZLACPY( 'Upper', P, N, BF, LDB, T, LDB )
*
*     Compute R - A*Q'
*
      CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, N, N, -CONE,
     $            A, LDA, Q, LDA, CONE, R, LDA )
*
*     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = ZLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) /
     $                 ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Q - Z'*B
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', P, N, P, CONE,
     $            Z, LDB, B, LDB, CZERO, BWK, LDB )
      CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, CONE, T, LDB,
     $            Q, LDA, -CONE, BWK, LDB )
*
*     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = ZLANGE( '1', P, N, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, M ) ) ) / BNORM ) /
     $                 ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL ZLASET( 'Full', N, N, CZERO, CONE, R, LDA )
      CALL ZHERK( 'Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R,
     $            LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL ZLASET( 'Full', P, P, CZERO, CONE, T, LDB )
      CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB,
     $            ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
      RETURN
*
*     End of ZGRQTS
*
      END
