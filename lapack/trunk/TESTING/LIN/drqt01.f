      SUBROUTINE DRQT01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
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
*  DRQT01 tests DGERQF, which computes the RQ factorization of an m-by-n
*  matrix A, and partially tests DORGRQ which forms the n-by-n
*  orthogonal matrix Q.
*
*  DRQT01 compares R with A*Q', and checks that Q is orthogonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the RQ factorization of A, as returned by DGERQF.
*          See DGERQF for further details.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The n-by-n orthogonal matrix Q.
*
*  R       (workspace) DOUBLE PRECISION array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L.
*          LDA >= max(M,N).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by DGERQF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(M,N))
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
      INTEGER            INFO, MINMN
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGERQF, DLACPY, DLASET, DORGRQ, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'DGERQF'
      CALL DGERQF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N )
     $      CALL DLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )
         IF( M.GT.1 )
     $      CALL DLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA,
     $                   Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 )
     $      CALL DLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA,
     $                   Q( 2, 1 ), LDA )
      END IF
*
*     Generate the n-by-n matrix Q
*
      SRNAMT = 'DORGRQ'
      CALL DORGRQ( N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL DLASET( 'Full', M, N, ZERO, ZERO, R, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 )
     $      CALL DLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA,
     $                   R( 1, N-M+1 ), LDA )
      ELSE
         IF( M.GT.N .AND. N.GT.0 )
     $      CALL DLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         IF( N.GT.0 )
     $      CALL DLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA,
     $                   R( M-N+1, 1 ), LDA )
      END IF
*
*     Compute R - A*Q'
*
      CALL DGEMM( 'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q,
     $            LDA, ONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
      RESID = DLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of DRQT01
*
      END
