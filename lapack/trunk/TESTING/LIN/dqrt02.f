      SUBROUTINE DQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
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
*  DQRT02 tests DORGQR, which generates an m-by-n matrix Q with
*  orthonornmal columns that is defined as the product of k elementary
*  reflectors.
*
*  Given the QR factorization of an m-by-n matrix A, DQRT02 generates
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A which was factorized by DQRT01.
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the QR factorization of A, as returned by DGEQRF.
*          See DGEQRF for further details.
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (LDA,N)
*
*  R       (workspace) DOUBLE PRECISION array, dimension (LDA,N)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R. LDA >= M.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (N)
*          The scalar factors of the elementary reflectors corresponding
*          to the QR factorization in AF.
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
*          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      EXTERNAL           DGEMM, DLACPY, DLASET, DORGQR, DSYRK
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
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the first k columns of the factorization to the array Q
*
      CALL DLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL DLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the first n columns of the matrix Q
*
      SRNAMT = 'DORGQR'
      CALL DORGQR( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R(1:n,1:k)
*
      CALL DLASET( 'Full', N, K, ZERO, ZERO, R, LDA )
      CALL DLACPY( 'Upper', N, K, AF, LDA, R, LDA )
*
*     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
*
      CALL DGEMM( 'Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A,
     $            LDA, ONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', M, K, A, LDA, RWORK )
      RESID = DLANGE( '1', N, K, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, R,
     $            LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of DQRT02
*
      END
