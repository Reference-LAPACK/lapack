      SUBROUTINE DRQT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
     $                   R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DRQT02 tests DORGRQ, which generates an m-by-n matrix Q with
*  orthonornmal rows that is defined as the product of k elementary
*  reflectors.
*
*  Given the RQ factorization of an m-by-n matrix A, DRQT02 generates
*  the orthogonal matrix Q defined by the factorization of the last k
*  rows of A; it compares R(m-k+1:m,n-m+1:n) with
*  A(m-k+1:m,1:n)*Q(n-m+1:n,1:n)', and checks that the rows of Q are
*  orthonormal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q to be generated.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q to be generated.
*          N >= M >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A which was factorized by DRQT01.
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the RQ factorization of A, as returned by DGERQF.
*          See DGERQF for further details.
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (LDA,N)
*
*  R       (workspace) DOUBLE PRECISION array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L. LDA >= N.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (M)
*          The scalar factors of the elementary reflectors corresponding
*          to the RQ factorization in AF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( R - A*Q' ) / ( N * norm(A) * EPS )
*          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ROGUE
      PARAMETER          ( ROGUE = -1.0D+10 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, DORGRQ, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the last k rows of the factorization to the array Q
*
      CALL DLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      IF( K.LT.N )
     $   CALL DLACPY( 'Full', K, N-K, AF( M-K+1, 1 ), LDA,
     $                Q( M-K+1, 1 ), LDA )
      IF( K.GT.1 )
     $   CALL DLACPY( 'Lower', K-1, K-1, AF( M-K+2, N-K+1 ), LDA,
     $                Q( M-K+2, N-K+1 ), LDA )
*
*     Generate the last n rows of the matrix Q
*
      SRNAMT = 'DORGRQ'
      CALL DORGRQ( M, N, K, Q, LDA, TAU( M-K+1 ), WORK, LWORK, INFO )
*
*     Copy R(m-k+1:m,n-m+1:n)
*
      CALL DLASET( 'Full', K, M, ZERO, ZERO, R( M-K+1, N-M+1 ), LDA )
      CALL DLACPY( 'Upper', K, K, AF( M-K+1, N-K+1 ), LDA,
     $             R( M-K+1, N-K+1 ), LDA )
*
*     Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'
*
      CALL DGEMM( 'No transpose', 'Transpose', K, M, N, -ONE,
     $            A( M-K+1, 1 ), LDA, Q, LDA, ONE, R( M-K+1, N-M+1 ),
     $            LDA )
*
*     Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', K, N, A( M-K+1, 1 ), LDA, RWORK )
      RESID = DLANGE( '1', K, M, R( M-K+1, N-M+1 ), LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL DLASET( 'Full', M, M, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', M, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of DRQT02
*
      END
