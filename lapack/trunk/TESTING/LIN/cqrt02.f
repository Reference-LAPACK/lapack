      SUBROUTINE CQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
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
      COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
     $                   R( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CQRT02 tests CUNGQR, which generates an m-by-n matrix Q with
*  orthonornmal columns that is defined as the product of k elementary
*  reflectors.
*
*  Given the QR factorization of an m-by-n matrix A, CQRT02 generates
*  the orthogonal matrix Q defined by the factorization of the first k
*  columns of A; it compares R(1:n,1:k) with Q(1:m,1:n)'*A(1:m,1:k),
*  and checks that the columns of Q are orthonormal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q to be generated.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q to be generated.
*          M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The m-by-n matrix A which was factorized by CQRT01.
*
*  AF      (input) COMPLEX array, dimension (LDA,N)
*          Details of the QR factorization of A, as returned by CGEQRF.
*          See CGEQRF for further details.
*
*  Q       (workspace) COMPLEX array, dimension (LDA,N)
*
*  R       (workspace) COMPLEX array, dimension (LDA,N)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R. LDA >= M.
*
*  TAU     (input) COMPLEX array, dimension (N)
*          The scalar factors of the elementary reflectors corresponding
*          to the QR factorization in AF.
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
*          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      EXTERNAL           CGEMM, CHERK, CLACPY, CLASET, CUNGQR
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
*     Copy the first k columns of the factorization to the array Q
*
      CALL CLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL CLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the first n columns of the matrix Q
*
      SRNAMT = 'CUNGQR'
      CALL CUNGQR( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R(1:n,1:k)
*
      CALL CLASET( 'Full', N, K, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA )
      CALL CLACPY( 'Upper', N, K, AF, LDA, R, LDA )
*
*     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', N, K, M,
     $            CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', M, K, A, LDA, RWORK )
      RESID = CLANGE( '1', N, K, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA,
     $            ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of CQRT02
*
      END
