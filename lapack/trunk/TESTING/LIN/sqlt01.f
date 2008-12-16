      SUBROUTINE SQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
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
      REAL               A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SQLT01 tests SGEQLF, which computes the QL factorization of an m-by-n
*  matrix A, and partially tests SORGQL which forms the m-by-m
*  orthogonal matrix Q.
*
*  SQLT01 compares L with Q'*A, and checks that Q is orthogonal.
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
*  A       (input) REAL array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) REAL array, dimension (LDA,N)
*          Details of the QL factorization of A, as returned by SGEQLF.
*          See SGEQLF for further details.
*
*  Q       (output) REAL array, dimension (LDA,M)
*          The m-by-m orthogonal matrix Q.
*
*  L       (workspace) REAL array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  TAU     (output) REAL array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by SGEQLF.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( L - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               ROGUE
      PARAMETER          ( ROGUE = -1.0E+10 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, MINMN
      REAL               ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SGEQLF, SLACPY, SLASET, SORGQL, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
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
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'SGEQLF'
      CALL SGEQLF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL SLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      IF( M.GE.N ) THEN
         IF( N.LT.M .AND. N.GT.0 )
     $      CALL SLACPY( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA )
         IF( N.GT.1 )
     $      CALL SLACPY( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA,
     $                   Q( M-N+1, M-N+2 ), LDA )
      ELSE
         IF( M.GT.1 )
     $      CALL SLACPY( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA,
     $                   Q( 1, 2 ), LDA )
      END IF
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'SORGQL'
      CALL SORGQL( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy L
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, L, LDA )
      IF( M.GE.N ) THEN
         IF( N.GT.0 )
     $      CALL SLACPY( 'Lower', N, N, AF( M-N+1, 1 ), LDA,
     $                   L( M-N+1, 1 ), LDA )
      ELSE
         IF( N.GT.M .AND. M.GT.0 )
     $      CALL SLACPY( 'Full', M, N-M, AF, LDA, L, LDA )
         IF( M.GT.0 )
     $      CALL SLACPY( 'Lower', M, M, AF( 1, N-M+1 ), LDA,
     $                   L( 1, N-M+1 ), LDA )
      END IF
*
*     Compute L - Q'*A
*
      CALL SGEMM( 'Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A,
     $            LDA, ONE, L, LDA )
*
*     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = SLANGE( '1', M, N, A, LDA, RWORK )
      RESID = SLANGE( '1', M, N, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, L, LDA )
      CALL SSYRK( 'Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, L,
     $            LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = SLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of SQLT01
*
      END
