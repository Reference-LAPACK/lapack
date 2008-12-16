      SUBROUTINE CLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK,
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
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CLQT02 tests CUNGLQ, which generates an m-by-n matrix Q with
*  orthonornmal rows that is defined as the product of k elementary
*  reflectors.
*
*  Given the LQ factorization of an m-by-n matrix A, CLQT02 generates
*  the orthogonal matrix Q defined by the factorization of the first k
*  rows of A; it compares L(1:k,1:m) with A(1:k,1:n)*Q(1:m,1:n)', and
*  checks that the rows of Q are orthonormal.
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
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The m-by-n matrix A which was factorized by CLQT01.
*
*  AF      (input) COMPLEX array, dimension (LDA,N)
*          Details of the LQ factorization of A, as returned by CGELQF.
*          See CGELQF for further details.
*
*  Q       (workspace) COMPLEX array, dimension (LDA,N)
*
*  L       (workspace) COMPLEX array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L. LDA >= N.
*
*  TAU     (input) COMPLEX array, dimension (M)
*          The scalar factors of the elementary reflectors corresponding
*          to the LQ factorization in AF.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( L - A*Q' ) / ( N * norm(A) * EPS )
*          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      REAL               ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      EXTERNAL           CLANGE, CLANSY, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CHERK, CLACPY, CLASET, CUNGLQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, REAL
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the first k rows of the factorization to the array Q
*
      CALL CLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL CLACPY( 'Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )
*
*     Generate the first n columns of the matrix Q
*
      SRNAMT = 'CUNGLQ'
      CALL CUNGLQ( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy L(1:k,1:m)
*
      CALL CLASET( 'Full', K, M, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA )
      CALL CLACPY( 'Lower', K, M, AF, LDA, L, LDA )
*
*     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', K, M, N,
     $            CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA )
*
*     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', K, N, A, LDA, RWORK )
      RESID = CLANGE( '1', K, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), L, LDA )
      CALL CHERK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of CLQT02
*
      END
